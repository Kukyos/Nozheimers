"""Build the results site from the pipeline outputs.

Run after pipeline.py:  python site/build.py
Writes site/data.js and copies structures and figures into site/. Every number on the page comes from here.
"""
import json, shutil
from pathlib import Path

import pandas as pd
import requests
from Bio.SeqUtils import seq1

ROOT = Path(__file__).resolve().parent.parent
SITE = ROOT / 'site'
GENES = {'APP': 351, 'PSEN1': 5663, 'PSEN2': 5664}
SPECIES = [('Human', 'Homo sapiens', 9606), ('Mouse', 'Mus musculus', 10090), ('Rat', 'Rattus norvegicus', 10116),
           ('Dog', 'Canis lupus familiaris', 9615), ('Cattle', 'Bos taurus', 9913)]
# experimental structures: chain -> (gene, offset added to the PDB residue number to get the human reference position)
STRUCTS = [
    dict(id='6IYC', label='PSEN1 gripping APP', chains={'B': ('PSEN1', 0), 'E': ('APP', 687)},
         note='Engineered changes for crosslinking: PSEN1 Q112C and APP V695C (SEQADV records).'),
    dict(id='7Y5Z', label='PSEN2 complex', chains={'B': ('PSEN2', 0)}, note=''),
    dict(id='5OQV', label='Amyloid-beta 42 fibril', chains={c: ('APP', 671) for c in 'ABCDEFGHI'},
         note='Nine stacked copies of amyloid-beta 1-42, the building block of plaques.'),
]
TOPO = {'Transmembrane', 'Intramembrane', 'Topological domain', 'Active site', 'Domain'}


def read(path, **kw):
    p = ROOT / path
    return pd.read_csv(p, **kw) if p.exists() else None


def ca_residues(pdb_text, chain):
    out = {}
    for l in pdb_text.splitlines():
        if l.startswith('ATOM') and l[21] == chain and l[12:16].strip() == 'CA':
            out[int(l[22:26])] = seq1(l[17:20])
    return out


def main():
    data = dict(genes={}, species=[dict(common=c, latin=l, taxid=t) for c, l, t in SPECIES], structures=[], figures=[])
    meta = pd.read_csv(ROOT / '01_metadata/metadata.csv')
    inc = meta[meta.decision.str.startswith('included')]
    data['sequences'] = inc[['gene', 'species', 'common_name', 'accession_version', 'transcript_protein_id', 'length',
                             'database_url', 'date_accessed', 'notes']].to_dict('records')
    data['n_candidates'] = int((meta.sequence_type == 'protein').sum())

    ident = read('05_conservation/identity_table.csv')
    blast = read('03_blast/blast_qc_summary.csv')
    allv = read('07_variants/clinvar_all_missense.csv')
    vmap = read('07_variants/variant_mapping.csv')
    data['identity'] = ident.to_dict('records') if ident is not None else []
    data['blast'] = blast.fillna('').to_dict('records') if blast is not None else []
    data['variants'] = vmap.fillna('').to_dict('records') if vmap is not None else []

    for gene, gid in GENES.items():
        g = dict(gene_id=gid)
        cons = read(f'05_conservation/{gene}_per_position_conservation.csv')
        if cons is not None:
            g['length'] = len(cons)
            g['residues'] = {c: ''.join(cons[c].astype(str)) for c, _, _ in SPECIES if c in cons}
            g['pct_fully_conserved'] = round(100 * cons.fully_conserved.mean(), 1)
        feats = next(iter((ROOT / '08_structures').glob(f'{gene}_*_uniprot_features.csv')), None)
        if feats is not None:
            g['uniprot'] = feats.stem.split('_')[1]
            f = pd.read_csv(feats)
            g['features'] = f[f.type.isin(TOPO)].fillna('').to_dict('records')
        if allv is not None:
            a = allv[allv.gene == gene]
            g['clinvar_total'] = len(a)
            g['clinvar_by_category'] = a.category.value_counts().to_dict()
            g['clinvar_selected'] = int((a.selection == 'selected').sum())
        if vmap is not None and cons is not None:
            m = vmap[(vmap.gene == gene) & (vmap.mapping == 'mapped')]
            pos = m.residue_position.unique()
            g['variant_positions'] = int(len(pos))
            g['variant_positions_conserved'] = int(cons[cons.human_position.isin(pos)].fully_conserved.sum())
        data['genes'][gene] = g

    pdb_dir = SITE / 'structures'
    pdb_dir.mkdir(exist_ok=True)
    for s in STRUCTS:
        src = ROOT / f'08_structures/pdb/{s["id"]}.pdb'
        text = src.read_text()
        (pdb_dir / src.name).write_text('\n'.join(l for l in text.splitlines() if l.startswith(('ATOM', 'TER', 'END'))) + '\n')
        entry_path = ROOT / f'08_structures/pdb/{s["id"]}_entry.json'
        if not entry_path.exists():
            entry_path.write_text(requests.get(f'https://data.rcsb.org/rest/v1/core/entry/{s["id"]}', timeout=60).text)
        e = json.loads(entry_path.read_text())
        resolved = {}
        for ch, (gene, off) in s['chains'].items():
            resolved.setdefault(gene, set()).update(k + off for k in ca_residues(text, ch))
        data['structures'].append(dict(
            s, file=f'structures/{src.name}', title=e['struct']['title'], method=e['exptl'][0]['method'].title(),
            resolution=(e.get('rcsb_entry_info', {}).get('resolution_combined') or [None])[0],
            released=e['rcsb_accession_info']['initial_release_date'][:10],
            resolved={g: sorted(v) for g, v in resolved.items()}))

    for af in sorted((ROOT / '08_structures').glob('AF-*.pdb')):
        acc = af.stem.split('-')[1]
        gene = next((k for k, v in data['genes'].items() if v.get('uniprot') == acc), None)
        if not gene:
            continue
        shutil.copy(af, pdb_dir / af.name)
        n = len(ca_residues(af.read_text(), 'A'))
        data['structures'].append(dict(id=af.stem, label=f'{gene} full length (AlphaFold prediction)', chains={'A': (gene, 0)},
                                       file=f'structures/{af.name}', title=f'AlphaFold DB model of human {gene} ({acc})',
                                       method='Predicted (AlphaFold)', resolution=None, released='',
                                       resolved={gene: list(range(1, n + 1))},
                                       note='A computer prediction, not an experiment. Colour by confidence is not shown here.'))

    fig_dir = SITE / 'figures'
    fig_dir.mkdir(exist_ok=True)
    for png in sorted((ROOT / '09_figures').glob('*.png')):
        shutil.copy(png, fig_dir / png.name)
        data['figures'].append(f'figures/{png.name}')

    steps = [('Download', ROOT / '01_metadata/metadata.csv'), ('BLAST check', ROOT / '03_blast/blast_qc_summary.csv'),
             ('Align', ROOT / '04_alignments'), ('Measure', ROOT / '05_conservation/identity_table.csv'),
             ('Trees', ROOT / '06_phylogeny'), ('Mutations', ROOT / '07_variants/variant_mapping.csv'),
             ('Structures', ROOT / '08_structures'), ('Report', ROOT / '10_report/report.docx')]
    data['status'] = [dict(step=n, done=p.exists() and (p.is_file() or any(p.iterdir()))) for n, p in steps]
    data['blast_pending'] = sorted({r['gene'] for r in data['blast'] if str(r.get('top_hit', '')).startswith('pending')})
    data['built'] = pd.Timestamp.now().strftime('%Y-%m-%d')
    (SITE / 'data.js').write_text('window.DATA = ' + json.dumps(data, separators=(',', ':'), default=str) + ';\n', encoding='utf-8')
    print('site/data.js written:', ', '.join(f'{k}={len(v) if hasattr(v, "__len__") else v}' for k, v in data.items()))


if __name__ == '__main__':
    main()
