"""Precompute the model-organism explorer (site/explorer/).

Run after pipeline.py:  python site/build_explorer.py
For each gene and each preset lab organism it finds the ortholog, aligns it to the human reference, superposes its
AlphaFold model on the human AlphaFold model, and scores how well the organism keeps the human protein: disease
residues, catalytic/amyloid-beta core, and residues that touch drugs in experimental human structures.
Everything is cached under 08_structures/ and 01_metadata/api_raw/ with the same rules as pipeline.py.
"""
import json, time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.SeqUtils import seq1
from Bio.SVDSuperimposer import SVDSuperimposer

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'site' / 'explorer'
API, MODELS, PDB = ROOT / '01_metadata/api_raw', ROOT / '08_structures/models', ROOT / '08_structures/pdb'
H = {'User-Agent': 'nozheimers_pipeline (Mozilla/5.0)'}
GENES = {'APP': (351, 'P05067'), 'PSEN1': (5663, 'P49768'), 'PSEN2': (5664, 'P49810')}
PRESETS = [  # taxid, common name, group
    (9598, 'Chimpanzee', 'Primates'), (9544, 'Rhesus macaque', 'Primates'),
    (10090, 'Mouse', 'Rodents'), (10116, 'Rat', 'Rodents'),
    (9615, 'Dog', 'Other mammals'), (9823, 'Pig', 'Other mammals'), (9913, 'Cattle', 'Other mammals'),
    (9031, 'Chicken', 'Birds'), (8364, 'Western clawed frog', 'Amphibians'), (7955, 'Zebrafish', 'Fish'),
    (7227, 'Fruit fly', 'Invertebrates'), (6239, 'C. elegans', 'Invertebrates'),
]
# NCBI makes no ortholog call for these; reviewed UniProt homologs (single presenilin = co-ortholog of PSEN1 and PSEN2)
HOMOLOGS = {7227: {'APP': ('P14599', 'Appl'), 'PSEN1': ('O02194', 'Psn'), 'PSEN2': ('O02194', 'Psn')},
            6239: {'APP': ('Q10651', 'apl-1'), 'PSEN1': ('P52166', 'sel-12'), 'PSEN2': ('P52166', 'sel-12')}}
# experimental human structures with a drug bound: (pdb, ligand code or partner chains, drug name, gene, kind)
DRUGS = [('6LQG', 'EN9', 'Avagacestat', 'PSEN1', 'gamma-secretase inhibitor'),
         ('6LR4', 'ESF', 'Semagacestat', 'PSEN1', 'gamma-secretase inhibitor'),
         ('7C9I', 'FTO', 'L-685,458', 'PSEN1', 'transition-state inhibitor'),
         ('7D8X', 'GZR', 'E2012', 'PSEN1', 'gamma-secretase modulator'),
         ('7Y5T', 'IGD', 'MRK-560', 'PSEN1', 'PSEN1-selective inhibitor'),
         ('6CO3', 'HL', 'Aducanumab', 'APP', 'anti-amyloid antibody'),
         ('4XXD', 'AB', 'Solanezumab', 'APP', 'anti-amyloid antibody'),
         ('5VZY', 'HL', 'Crenezumab', 'APP', 'anti-amyloid antibody'),
         ('5CSZ', 'HL', 'Gantenerumab', 'APP', 'anti-amyloid antibody')]
CONTACT = 4.5  # Å, heavy atoms
CONF = 70  # pLDDT cut-off for residues used in the superposition
AL = PairwiseAligner(mode='global', substitution_matrix=substitution_matrices.load('BLOSUM62'),
                     open_gap_score=-11, extend_gap_score=-1)
B62 = substitution_matrices.load('BLOSUM62')


def get(url, **kw):
    for i in range(6):
        try:
            r = requests.get(url, headers=H, timeout=120, **kw)
            r.raise_for_status()
            return r.text
        except requests.RequestException:
            if i == 5:
                raise
            time.sleep(5 * (i + 1))


def cached(path, url, **kw):
    path = Path(path)
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(get(url, **kw), encoding='utf-8')
    return path.read_text(encoding='utf-8')


def uniprot(acc):
    return json.loads(cached(API / f'uniprot_{acc}.json', f'https://rest.uniprot.org/uniprotkb/{acc}.json'))


def alphafold(acc):
    meta = json.loads(cached(MODELS / f'alphafold_{acc}_prediction.json', f'https://alphafold.ebi.ac.uk/api/prediction/{acc}'))[0]
    pdb = cached(MODELS / Path(meta['pdbUrl']).name, meta['pdbUrl'])
    return meta, pdb


def atoms(pdb_text, chain=None):
    """[(chain, resi, resn, atom name, xyz, bfactor)] for ATOM/HETATM records."""
    out = []
    for l in pdb_text.splitlines():
        if l.startswith(('ATOM', 'HETATM')) and (chain is None or l[21] in chain):
            out.append((l[21], int(l[22:26]), l[17:20].strip(), l[12:16].strip(), np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]), float(l[60:66]), l[:6].strip()))
    return out


def ca_map(pdb_text):
    return {a[1]: (a[4], a[5]) for a in atoms(pdb_text) if a[3] == 'CA' and a[6] == 'ATOM'}


def align(h, o):
    """human position -> organism position (1-based), for aligned residue pairs."""
    a = AL.align(h, o)[0]
    m = {}
    for (hs, he), (os_, oe) in zip(*a.aligned):
        for k in range(he - hs):
            m[int(hs + k + 1)] = int(os_ + k + 1)
    return m, a.score


def ortholog_reports(gid):
    return {int(r['gene']['tax_id']): r['gene'] for r in json.loads((API / f'gene_{gid}_orthologs_all.json').read_text())['reports']}


def pick_ortholog(gene, tax, human_seq, reports):
    if tax in HOMOLOGS:
        acc, sym = HOMOLOGS[tax][gene]
        note = 'No NCBI ortholog call; reviewed UniProt homolog.'
        if gene.startswith('PSEN'):
            note += ' This organism has one presenilin, a co-ortholog of both PSEN1 and PSEN2.'
        return dict(acc=acc, symbol=sym, gene_id=None, source='UniProt homolog', note=note)
    g = reports.get(tax)
    if not g:
        return None
    # every UniProt entry linked to the NCBI GeneID (Swiss-Prot entries can be fragments); best alignment score wins
    res = json.loads(cached(API / f'uniprot_geneid_{g["gene_id"]}.json', 'https://rest.uniprot.org/uniprotkb/search',
                            params=dict(query=f'xref:geneid-{g["gene_id"]}', fields='accession', format='json', size=25)))['results']
    accs = list(dict.fromkeys(list(g.get('swiss_prot_accessions') or []) + [r['primaryAccession'] for r in res]))
    if not accs:
        return None
    scored = sorted(accs, key=lambda a: -AL.score(human_seq, uniprot(a)['sequence']['value']))
    return dict(acc=scored[0], symbol=g['symbol'], gene_id=int(g['gene_id']), source='NCBI ortholog',
                note='' if len(accs) == 1 else f'{len(accs)} UniProt entries for this gene; best alignment to human kept.')


def drug_contacts():
    """human residue positions within CONTACT Å of each drug, from the experimental complex."""
    refs = {g: uniprot(acc)['sequence']['value'] for g, (_, acc) in GENES.items()}
    out = []
    for pdb_id, lig, name, gene, kind in DRUGS:
        text = (PDB / f'{pdb_id}.pdb').read_text()
        at = [a for a in atoms(text) if not a[3].startswith('H')]
        if gene == 'PSEN1':
            target = [a for a in at if a[0] == 'B' and a[6] == 'ATOM']
            partner = [a for a in at if a[2] == lig]
            offset = 0
        else:  # amyloid-beta peptide chain(s) vs antibody chains; find the offset by matching the peptide to human APP
            partner = [a for a in at if a[0] in lig and a[6] == 'ATOM']
            pep_chains = sorted({a[0] for a in at if a[0] not in lig and a[6] == 'ATOM'})
            pep_chain = next(c for c in pep_chains if len({a[1] for a in at if a[0] == c}) < 60)
            target = [a for a in at if a[0] == pep_chain and a[6] == 'ATOM']
            res = sorted({(a[1], a[2]) for a in target})
            pep = ''.join(seq1(n) for _, n in res)
            i = refs['APP'].find(pep)
            assert i >= 0, f'{pdb_id}: peptide {pep} not found in human APP'
            offset = i + 1 - res[0][0]
        assert partner, f'{pdb_id}: no atoms for {lig}'
        P = np.array([a[4] for a in partner])
        hits = set()
        for a in target:
            if np.min(np.linalg.norm(P - a[4], axis=1)) <= CONTACT:
                hits.add(a[1] + offset)
        hits = sorted(hits)
        mism = [f'{refs[gene][p - 1]}{p}{seq1(next(a[2] for a in target if a[1] + offset == p))}' for p in hits
                if refs[gene][p - 1] != seq1(next(a[2] for a in target if a[1] + offset == p))]
        assert len(mism) <= 1, f'{pdb_id}: numbering does not match human {gene}: {mism}'  # allows one engineered residue
        meta = json.loads(cached(PDB / f'{pdb_id}_entry.json', f'https://data.rcsb.org/rest/v1/core/entry/{pdb_id}'))
        out.append(dict(pdb=pdb_id, ligand=lig, drug=name, gene=gene, kind=kind, residues=hits,
                        human=''.join(refs[gene][p - 1] for p in hits),
                        method=meta['exptl'][0]['method'].title(),
                        resolution=(meta.get('rcsb_entry_info', {}).get('resolution_combined') or [None])[0]))
        print(f'  {pdb_id} {name}: {len(hits)} {gene} residues within {CONTACT} A')
    return out


def superpose(human_ca, org_ca, m, h_pl, o_pl, regions):
    """Whole-protein fit if it is tight; otherwise fit on the best-fitting domain (multi-domain proteins like APP
    have flexible linkers, so one rigid fit over all domains is meaningless)."""
    fits = [dict(kabsch(human_ca, org_ca, m, h_pl, o_pl) or {}, region='whole protein')]
    for label, lo, hi in regions:
        f = kabsch(human_ca, org_ca, m, h_pl, o_pl, (lo, hi))
        if f:
            fits.append(dict(f, region=f'{label} ({lo}-{hi})'))
    fits = [f for f in fits if f.get('n')]
    if not fits:
        return None
    best = fits[0] if fits[0]['region'] == 'whole protein' and fits[0]['rmsd'] <= 3 else \
        max([f for f in fits if f['rmsd'] <= 3] or fits, key=lambda f: (f['rmsd'] <= 3, f['n'], -f['rmsd']))
    return dict(best, fits=[dict(region=f['region'], rmsd=f['rmsd'], n=f['n']) for f in fits])


def kabsch(human_ca, org_ca, m, h_pl, o_pl, span=None):
    pairs = [(hp, op) for hp, op in m.items() if hp in human_ca and op in org_ca
             and h_pl.get(hp, 0) >= CONF and o_pl.get(op, 0) >= CONF and (span is None or span[0] <= hp <= span[1])]
    if len(pairs) < 10:
        return None
    X = np.array([human_ca[hp][0] for hp, _ in pairs])
    Y = np.array([org_ca[op][0] for _, op in pairs])
    sup = SVDSuperimposer()
    sup.set(X, Y)
    sup.run()
    rot, tran = sup.get_rotran()
    return dict(rmsd=round(float(sup.get_rms()), 2), n=len(pairs), rot=rot.round(6).tolist(), tran=tran.round(4).tolist())


GROUPS = [('Primates', 'Primates'), ('Rodentia', 'Rodents'), ('Lagomorpha', 'Other mammals'), ('Mammalia', 'Other mammals'),
          ('Aves', 'Birds'), ('Lepidosauria', 'Reptiles'), ('Testudines', 'Reptiles'), ('Crocodylia', 'Reptiles'),
          ('Amphibia', 'Amphibians'), ('Actinopterygii', 'Fish'), ('Chondrichthyes', 'Fish'), ('Sarcopterygii', 'Fish'),
          ('Cyclostomata', 'Fish')]


def write_orthologs():
    """Every organism NCBI lists an ortholog for, per gene, with an organism type from the UniProt taxonomy lineage."""
    reports = {g: ortholog_reports(gid) for g, (gid, _) in GENES.items()}
    taxa = sorted({t for r in reports.values() for t in r})
    lin_path = API / 'uniprot_taxonomy_lineages.json'
    lin = json.loads(lin_path.read_text()) if lin_path.exists() else {}
    todo = [t for t in taxa if str(t) not in lin]
    for i in range(0, len(todo), 100):
        q = ' OR '.join(f'tax_id:{t}' for t in todo[i:i + 100])
        res = json.loads(get('https://rest.uniprot.org/taxonomy/search', params=dict(query=q, fields='id,lineage', format='json', size=500)))['results']
        for r in res:
            lin[str(r['taxonId'])] = [x['scientificName'] for x in r.get('lineage', [])]
        lin_path.write_text(json.dumps(lin))
    def group(t):
        names = set(lin.get(str(t), []))
        return next((label for clade, label in GROUPS if clade in names), 'Other')
    out = {g: [[t, r['taxname'], r.get('common_name', ''), int(r['gene_id']), r['symbol'], r.get('swiss_prot_accessions') or [], group(t)]
               for t, r in sorted(rep.items(), key=lambda x: x[1]['taxname'])] for g, rep in reports.items()}
    (OUT / 'orthologs.js').write_text('window.ORTHOLOGS = ' + json.dumps(out, separators=(',', ':')) + ';' + chr(10), encoding='utf-8')
    from collections import Counter
    print('orthologs:', {g: len(v) for g, v in out.items()}, Counter(x[6] for x in out['PSEN1']))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'models').mkdir(exist_ok=True)
    vm = pd.read_csv(ROOT / '07_variants/variant_mapping.csv')
    vm = vm[vm.mapping == 'mapped']
    drugs = drug_contacts()
    result = dict(genes={}, drugs=drugs, presets=[dict(taxid=t, name=n, group=g) for t, n, g in PRESETS],
                  criteria=dict(contact_A=CONTACT, plddt_cutoff=CONF), built=time.strftime('%Y-%m-%d'))
    for gene, (gid, hacc) in GENES.items():
        print(gene)
        hu = uniprot(hacc)
        hseq = hu['sequence']['value']
        hmeta, hpdb = alphafold(hacc)
        (OUT / 'models' / f'human_{gene}.pdb').write_text('\n'.join(l for l in hpdb.splitlines() if l.startswith(('ATOM', 'TER', 'END'))) + '\n')
        hca = ca_map(hpdb)
        h_pl = {k: v[1] for k, v in hca.items()}
        feats = [dict(type=f['type'], description=f.get('description', ''), start=f['location']['start']['value'],
                      end=f['location']['end']['value']) for f in hu['features']
                 if f['type'] in ('Transmembrane', 'Intramembrane', 'Topological domain', 'Active site', 'Domain', 'Region', 'Chain', 'Peptide')]
        disease = vm[vm.gene == gene]
        dis_pos = sorted(disease.residue_position.unique().tolist())
        core = list(range(672, 714)) if gene == 'APP' else sorted({p for f in feats if f['type'] == 'Transmembrane' for p in range(f['start'], f['end'] + 1)})
        active = [f['start'] for f in feats if f['type'] == 'Active site']
        gdrugs = [d for d in drugs if d['gene'] == gene]
        regions = [(f['description'] or f['type'], f['start'], f['end']) for f in feats if f['type'] == 'Domain' and f['end'] - f['start'] >= 30]
        if gene == 'APP':
            regions.append(('amyloid-beta and membrane helix', 672, 723))
        # sanity: human against itself
        m, _ = align(hseq, hseq)
        assert kabsch(hca, hca, m, h_pl, h_pl)['rmsd'] < 1e-3
        orgs = {}
        reports = ortholog_reports(gid)
        for tax, name, group in PRESETS:
            o = pick_ortholog(gene, tax, hseq, reports)
            if not o:
                orgs[tax] = dict(name=name, group=group, missing=f'No {gene} ortholog listed by NCBI for this organism.')
                print(f'  {name}: no ortholog')
                continue
            ou = uniprot(o['acc'])
            oseq = ou['sequence']['value']
            try:
                ometa, opdb = alphafold(o['acc'])
            except requests.RequestException:
                ometa, opdb = None, None
            m, score = align(hseq, oseq)
            same = [hp for hp, op in m.items() if hseq[hp - 1] == oseq[op - 1]]
            similar = [hp for hp, op in m.items() if hseq[hp - 1] != oseq[op - 1] and B62[hseq[hp - 1]][oseq[op - 1]] > 0]
            res_at = lambda hp: oseq[m[hp] - 1] if hp in m else '-'
            rec = dict(name=name, group=group, latin=ou['organism']['scientificName'], uniprot=o['acc'], symbol=o['symbol'],
                       gene_id=o['gene_id'], source=o['source'], note=o['note'], length=len(oseq),
                       identity=round(100 * len(same) / len(m), 1), coverage=round(100 * len(m) / len(hseq), 1),
                       similar=round(100 * (len(same) + len(similar)) / len(m), 1),
                       map={hp: m[hp] for hp in m}, seq=oseq,
                       disease_total=len(dis_pos), disease_same=sum(res_at(p) == hseq[p - 1] for p in dis_pos),
                       core_total=len(core), core_same=sum(res_at(p) == hseq[p - 1] for p in core),
                       active=[dict(pos=p, human=hseq[p - 1], org=res_at(p)) for p in active],
                       drugs=[dict(pdb=d['pdb'], drug=d['drug'], total=len(d['residues']),
                                   same=sum(res_at(p) == hseq[p - 1] for p in d['residues']),
                                   diffs=[f'{hseq[p - 1]}{p}{res_at(p)}' for p in d['residues'] if res_at(p) != hseq[p - 1]]) for d in gdrugs])
            if ometa:
                if ometa.get('uniprotSequence', ometa.get('sequence')) != oseq:
                    rec['note'] += ' AlphaFold model sequence differs from the current UniProt sequence.'
                fname = f'{gene}_{tax}.pdb'
                (OUT / 'models' / fname).write_text('\n'.join(l for l in opdb.splitlines() if l.startswith(('ATOM', 'TER', 'END'))) + '\n')
                oca = ca_map(opdb)
                o_pl = {k: v[1] for k, v in oca.items()}
                rec.update(model=f'explorer/models/{fname}', alphafold=ometa.get('entryId'),
                           plddt_mean=round(float(np.mean(list(o_pl.values()))), 1),
                           superpose=superpose(hca, oca, m, h_pl, o_pl, regions))
            orgs[tax] = rec
            sp = rec.get('superpose') or {}
            print(f'  {name}: {o["acc"]} {o["symbol"]} id {rec["identity"]}% disease {rec["disease_same"]}/{rec["disease_total"]} '
                  f'core {rec["core_same"]}/{rec["core_total"]} rmsd {sp.get("rmsd")} on {sp.get("n")} CA')
        result['genes'][gene] = dict(
            uniprot=hacc, gene_id=gid, seq=hseq, length=len(hseq), features=feats, alphafold=hmeta.get('entryId'),
            plddt=[h_pl.get(i, 0) for i in range(1, len(hseq) + 1)], human_model=f'explorer/models/human_{gene}.pdb',
            core_label='amyloid-beta (APP 672-713)' if gene == 'APP' else 'membrane-spanning helices',
            disease=[dict(pos=int(r.residue_position), hgvs=r.hgvs_p, id=int(r.variation_id), url=r.record_url,
                          cls=r.classification) for r in disease.itertuples()],
            organisms=orgs)
    write_orthologs()
    (OUT / 'explorer.js').write_text('window.EXPLORER = ' + json.dumps(result, separators=(',', ':')) + ';\n', encoding='utf-8')
    print('written', OUT / 'explorer.js', round((OUT / 'explorer.js').stat().st_size / 1e6, 2), 'MB')


if __name__ == '__main__':
    main()
