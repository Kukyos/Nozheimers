"""Comparative sequence analysis of APP, PSEN1 and PSEN2 across five mammals.

Implements Steps 1-12 of Alzheimers_Gene_Sequence_Comparative_Analysis_Project_Guide.docx.
Run:  python pipeline.py
Every download is cached on disk (with its first access date), so re-running is
fast and never overwrites an original. Delete a cached file to fetch it again.
"""
import io, json, random, re, time
from collections import Counter
from datetime import date, datetime
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import requests
from Bio import AlignIO, Phylo, SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio.SeqUtils import seq1

EMAIL = 'amohamedarmaan@gmail.com'  # required by EBI Job Dispatcher, requested by NCBI
TOOL = 'nozheimers_pipeline'
GENES = {'APP': 351, 'PSEN1': 5663, 'PSEN2': 5664}  # NCBI Gene IDs from the guide
SPECIES = {9606: ('Homo sapiens', 'Human'), 10090: ('Mus musculus', 'Mouse'), 10116: ('Rattus norvegicus', 'Rat'),
           9615: ('Canis lupus familiaris', 'Dog'), 9913: ('Bos taurus', 'Cattle')}
BOOTSTRAPS, SEED = 1000, 1
NP_TOLERANCE = 0.95  # keep a curated NP_ isoform if its score is within 5% of the best candidate

ROOT = Path(__file__).resolve().parent
D = {n: ROOT / n for n in ['01_metadata', '02_raw_sequences', '03_blast', '04_alignments', '05_conservation',
                           '06_phylogeny', '07_variants', '08_structures', '09_figures', '10_report']}
API, ORIG = D['01_metadata'] / 'api_raw', D['02_raw_sequences'] / 'original'
EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
DATASETS = 'https://api.ncbi.nlm.nih.gov/datasets/v2'
BLAST_URL = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'
EBI = 'https://www.ebi.ac.uk/Tools/services/rest/clustalo'
NCBI = dict(tool=TOOL, email=EMAIL)
TODAY = date.today().isoformat()
AA = set('ACDEFGHIKLMNPQRSTVWY')
LOG = D['01_metadata'] / 'methods_log.md'
ACCESS_FILE = D['01_metadata'] / 'access_dates.json'
ACCESS = json.loads(ACCESS_FILE.read_text()) if ACCESS_FILE.exists() else {}
SETTINGS = []  # rows for the analysis-settings table


# ---------- helpers ----------
def log(msg):
    with LOG.open('a', encoding='utf-8') as f:
        f.write(f'- {msg}\n')
    print(msg)


def get(url, params=None, data=None, tries=6):
    for i in range(tries):
        try:
            r = requests.request('POST' if data else 'GET', url, params=params, data=data, timeout=180,
                                 headers={'User-Agent': f'{TOOL} ({EMAIL})'})
            r.raise_for_status()
            if 'ncbi.nlm.nih.gov' in url:
                time.sleep(0.4)  # ponytail: no API key, so stay under NCBI's 3 requests/s
            return r.text
        except requests.RequestException:
            if i == tries - 1:
                raise
            time.sleep(10 * (i + 1))


def cached(path, fetch):
    path = Path(path)
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        text = fetch()
        path.write_text(text, encoding='utf-8')
        ACCESS[path.relative_to(ROOT).as_posix()] = TODAY
        ACCESS_FILE.write_text(json.dumps(ACCESS, indent=1))
    return path.read_text(encoding='utf-8')


def accessed(path):
    return ACCESS.get(Path(path).relative_to(ROOT).as_posix(), TODAY)


def efetch(db, ids, rettype):
    return get(EUTILS + 'efetch.fcgi', data=dict(db=db, id=ids, rettype=rettype, retmode='text', **NCBI))


def records(fasta_text):
    """[(accession.version, header, sequence, original record text)] keeping the original text verbatim."""
    out = []
    for chunk in ('\n' + fasta_text.strip()).split('\n>')[1:]:
        header, *lines = chunk.split('\n')
        out.append((header.split()[0], header, ''.join(l.strip() for l in lines), '>' + chunk.strip() + '\n'))
    return out


def setting(step, tool, version, params):
    SETTINGS.append(dict(step=step, tool_or_resource=tool, version_or_endpoint=version, parameters=params,
                         date=TODAY))


def ident(a, b):
    """The project's single identity definition: identical / columns where both sequences have a residue."""
    both = [(x, y) for x, y in zip(a, b) if x != '-' and y != '-']
    return sum(x == y for x, y in both), len(both)


def pid(a, b):
    same, shared = ident(a, b)
    return same / shared if shared else 0.0


# ---------- Step 1: folders and log ----------
def setup():
    for p in list(D.values()) + [API, ORIG]:
        p.mkdir(parents=True, exist_ok=True)
    LOG.write_text(f'# Methods and decision log\n\nRun started {datetime.now():%Y-%m-%d %H:%M}. '
                   'Rules are fixed here before any retrieval; later entries record what was done.\n\n'
                   '## Predefined rules\n', encoding='utf-8')
    log('Genes: APP (Gene ID 351), PSEN1 (5663), PSEN2 (5664), each analysed separately. APOE not included '
        '(optional extension in the guide).')
    log('Species: ' + ', '.join(f'{s} ({c.lower()}, taxid {t})' for t, (s, c) in SPECIES.items()) +
        '. A species is excluded for a gene if NCBI lists no ortholog; exclusions are logged.')
    log('Human reference: the MANE Select transcript and its protein, read from the NCBI Gene product report.')
    log('Ortholog isoform rule: among all RefSeq proteins of the NCBI-designated ortholog gene, take the one with '
        'the highest global alignment score to the human reference (BLOSUM62, gap open -11, extend -1); a curated '
        f'NP_ record is preferred over a predicted XP_ record when its score is within {100 - NP_TOLERANCE * 100:.0f}% '
        'of the best.')
    log('Sequence type: protein for all alignment, identity, phylogeny and variant analyses. Human CDS is saved '
        'for reference only and not mixed into protein analyses.')
    log('Percent identity: identical residues / alignment columns where both sequences have a residue (gap '
        'columns excluded). Coverage: those shared columns / human protein length. The same definition gives '
        'the p-distance (1 - identity) used for the trees.')
    log(f'Trees: neighbour-joining on p-distance, {BOOTSTRAPS} bootstrap replicates (columns resampled, seed {SEED}), '
        'midpoint-rooted for display only.')
    log('ClinVar selection rule: missense single-gene records for the gene whose germline classification is '
        'Pathogenic, Likely pathogenic or Pathogenic/Likely pathogenic; review status at least "criteria provided, '
        'single submitter"; at least one condition name containing "Alzheimer". All other records are kept in '
        'the full table but not mapped.')
    log('Variant mapping rule: map only when the ClinVar HGVS transcript is the human reference transcript and the '
        'reference amino acid matches the reference protein at that position; otherwise label "unmapped".')
    with LOG.open('a', encoding='utf-8') as f:
        f.write('\n## Work log\n')


# ---------- Steps 2-4: sequences, metadata, validation ----------
def human_refs(meta):
    refs = {}
    for gene, gid in GENES.items():
        rep = json.loads(cached(API / f'gene_{gid}_products.json',
                                lambda: get(f'{DATASETS}/gene/id/{gid}/product_report')))['reports'][0]['product']
        mane = [t for t in rep['transcripts'] if t.get('select_category') == 'MANE_SELECT']
        assert len(mane) == 1, f'{gene}: expected exactly one MANE Select transcript'
        nm, npa = mane[0]['accession_version'], mane[0]['protein']['accession_version']
        pf, cf = ORIG / f'{gene}_Homo_sapiens_{npa}.fasta', ORIG / f'{gene}_Homo_sapiens_{nm}_CDS.fasta'
        (acc, header, seq, _), = records(cached(pf, lambda: efetch('protein', npa, 'fasta')))
        (_, cheader, cseq, _), = records(cached(cf, lambda: efetch('nuccore', nm, 'fasta_cds_na')))
        refs[gene] = dict(gid=gid, nm=nm, np=npa, seq=seq, header=header)
        log(f'{gene}: {len(rep["transcripts"])} RefSeq transcripts listed; selected MANE Select {nm} -> {npa} '
            f'({len(seq)} aa). Saved {pf.name} and {cf.name}.')
        common = dict(gene=gene, species='Homo sapiens', common_name='human', ncbi_gene_id=gid, date_accessed=accessed(pf),
                      decision='included (reference)')
        meta.append(dict(common, sequence_type='protein', accession_version=npa, transcript_protein_id=nm, length=len(seq),
                         database_url=f'https://www.ncbi.nlm.nih.gov/protein/{npa}', record_description=header,
                         notes='MANE Select protein', file=pf.name))
        meta.append(dict(common, sequence_type='CDS (nucleotide)', accession_version=nm, transcript_protein_id=npa,
                         length=len(cseq), database_url=f'https://www.ncbi.nlm.nih.gov/nuccore/{nm}',
                         record_description=cheader, notes='MANE Select CDS; saved for reference, not used in protein analyses',
                         decision='saved (not analysed)', file=cf.name))
    setting('2 Human references', 'NCBI Datasets v2 gene product_report; E-utilities efetch', DATASETS,
            'MANE Select transcript and protein per gene')
    return refs


def orthologs(refs, meta):
    aligner = PairwiseAligner(mode='global', substitution_matrix=substitution_matrices.load('BLOSUM62'),
                              open_gap_score=-11, extend_gap_score=-1)
    taxa = '&'.join(f'taxon_filter={t}' for t in SPECIES if t != 9606)
    seqs = {g: {9606: (r['np'], r['seq'])} for g, r in refs.items()}
    candidates = {g: {} for g in GENES}
    for gene, gid in GENES.items():
        orth = json.loads(cached(API / f'gene_{gid}_orthologs.json',
                                 lambda: get(f'{DATASETS}/gene/id/{gid}/orthologs?{taxa}')))
        by_tax = {int(r['gene']['tax_id']): r['gene'] for r in orth.get('reports', [])}
        for tax, (sci, common) in SPECIES.items():
            if tax == 9606:
                continue
            if tax not in by_tax:
                log(f'{gene}: EXCLUDED {sci}, no NCBI ortholog listed.')
                continue
            og = by_tax[tax]
            ogid = int(og['gene_id'])
            prod = json.loads(cached(API / f'gene_{ogid}_products.json',
                                     lambda: get(f'{DATASETS}/gene/id/{ogid}/product_report')))['reports'][0]['product']
            tx_of = {t['protein']['accession_version']: t['accession_version'] for t in prod['transcripts'] if t.get('protein')}
            sp = sci.replace(' ', '_')
            cf = D['02_raw_sequences'] / 'candidates' / f'{gene}_{sp}_candidate_proteins.fasta'
            recs = records(cached(cf, lambda: efetch('protein', ','.join(sorted(tx_of)), 'fasta')))
            scored = sorted(((aligner.score(refs[gene]['seq'], s), a, h, s, raw) for a, h, s, raw in recs),
                            key=lambda x: (-x[0], not x[1].startswith('NP_'), x[1]))
            best = scored[0]
            nps = [x for x in scored if x[1].startswith('NP_') and x[0] >= NP_TOLERANCE * best[0]]
            pick = nps[0] if nps else best
            score, acc, header, seq, raw = pick
            out = ORIG / f'{gene}_{sp}_{acc}.fasta'
            if not out.exists():
                out.write_text(raw, encoding='utf-8')
            seqs[gene][tax] = (acc, seq)
            candidates[gene][tax] = dict(gene_id=ogid, symbol=og['symbol'], proteins={a for a, *_ in recs})
            flags = [w for w, hit in [('predicted (XP_) record', acc.startswith('XP_')),
                                      ('LOW QUALITY PROTEIN in header', 'LOW QUALITY' in header),
                                      ('partial in header', 'partial' in header.lower())] if hit]
            log(f'{gene}: {sci} ortholog {og["symbol"]} (Gene ID {ogid}) has {len(recs)} RefSeq proteins; selected '
                f'{acc} ({len(seq)} aa, score {score:.0f}; best candidate score {best[0]:.0f})'
                + (f'; FLAGS: {", ".join(flags)}' if flags else '') + '.')
            for s, a, h, sq, _ in scored:
                chosen = a == acc
                meta.append(dict(gene=gene, species=sci, common_name=common.lower(), sequence_type='protein',
                                 ncbi_gene_id=ogid, accession_version=a, transcript_protein_id=tx_of.get(a, ''),
                                 length=len(sq), database_url=f'https://www.ncbi.nlm.nih.gov/protein/{a}',
                                 date_accessed=accessed(cf), record_description=h,
                                 decision='included' if chosen else 'excluded (other isoform)',
                                 notes=(f'alignment score to human {s:.0f}' + ('; ' + ', '.join(flags) if chosen and flags else '')
                                        + ('' if chosen else f'; selected isoform {acc} scores {score:.0f}')),
                                 file=out.name if chosen else cf.name))
    setting('3 Orthologs', 'NCBI Datasets v2 gene orthologs + product_report; Biopython PairwiseAligner',
            f'{DATASETS}/gene/id/{{id}}/orthologs', 'global alignment, BLOSUM62, gap open -11, extend -1; '
            f'NP_ preferred within {100 - NP_TOLERANCE * 100:.0f}% of best score')
    return seqs, candidates


def validate(seqs, meta):
    for gene, by_tax in seqs.items():
        accs = [a for a, _ in by_tax.values()]
        assert len(set(accs)) == len(accs), f'{gene}: duplicate accession'
        for tax, (acc, seq) in by_tax.items():
            bad = set(seq) - AA
            assert not bad - {'X', 'U'}, f'{acc}: unexpected characters {bad}'
            if bad:
                log(f'{gene}: {acc} contains {bad}; kept as downloaded and noted.')
            assert 0.5 < len(seq) / len(by_tax[9606][1]) < 2, f'{acc}: implausible length {len(seq)}'
    for row in meta:
        if row['sequence_type'] == 'protein':
            assert row['accession_version'] in row['record_description'], row
            assert f"[{row['species']}]" in row['record_description'], row
    log('FASTA validation passed: amino-acid alphabet only, no duplicate accessions, lengths within 0.5-2x of '
        'human, one sequence type (protein) per analysis set, headers carry accession and [species]. No cleaned '
        'copies were needed, so 02_raw_sequences/cleaned/ is not used.')
    cols = ['gene', 'species', 'common_name', 'sequence_type', 'ncbi_gene_id', 'accession_version',
            'transcript_protein_id', 'length', 'database_url', 'date_accessed', 'record_description', 'decision',
            'notes', 'file']
    df = pd.DataFrame(meta)[cols]
    df.to_csv(D['01_metadata'] / 'metadata.csv', index=False)
    df.to_excel(D['01_metadata'] / 'metadata.xlsx', index=False)
    return df


# ---------- Step 5: BLAST QC ----------
def blast_qc(refs, seqs, candidates):
    entrez = ' OR '.join(f'txid{t}[ORGN]' for t in SPECIES if t != 9606)
    rows = []
    for gene, ref in refs.items():
        path = D['03_blast'] / f'{gene}_blastp_{ref["np"]}_refseq_protein.json'

        def run():
            rid_file = path.with_suffix('.rid')  # resume a submitted search after a dropped connection
            if rid_file.exists():
                rid = rid_file.read_text().strip()
                log(f'{gene}: resuming BLASTp RID {rid}.')
            else:
                put = get(BLAST_URL, data=dict(CMD='Put', PROGRAM='blastp', DATABASE='refseq_protein', QUERY=ref['np'],
                                               ENTREZ_QUERY=entrez, HITLIST_SIZE=100, EMAIL=EMAIL, TOOL=TOOL))
                rid = re.search(r'RID = (\S+)', put).group(1)
                rid_file.write_text(rid)
                log(f'{gene}: BLASTp submitted, RID {rid}.')
            for _ in range(60):
                time.sleep(60)  # NCBI asks for at most one status poll per minute per RID
                try:
                    status = re.search(r'Status=(\w+)', get(BLAST_URL, params=dict(CMD='Get', FORMAT_OBJECT='SearchInfo', RID=rid))).group(1)
                except requests.RequestException:
                    continue  # ponytail: flaky network here; keep polling, the search keeps running at NCBI
                if status == 'READY':
                    text = get(BLAST_URL, params=dict(CMD='Get', FORMAT_TYPE='JSON2_S', RID=rid))
                    rid_file.unlink()
                    return text
                assert status == 'WAITING', f'BLAST {rid} status {status}'
            raise TimeoutError(f'BLAST {rid} not ready after 60 min; re-run to resume')

        report = json.loads(cached(path, run))['BlastOutput2'][0]['report']
        search = report['results']['search']
        for tax, (sci, common) in SPECIES.items():
            if tax == 9606 or tax not in seqs[gene]:
                continue
            hit = next((h for h in search['hits'] if any(d.get('taxid') == tax for d in h['description'])), None)
            if hit is None:
                rows.append(dict(gene=gene, species=sci, top_hit='none'))
                log(f'{gene}: WARNING no BLAST hit for {sci}.')
                continue
            hit_accs = {d['accession'] for d in hit['description'] if d.get('taxid') == tax}
            base = lambda a: a.split('.')[0]
            orth = {base(a) for a in candidates[gene][tax]['proteins']}
            hsp = hit['hsps'][0]
            covered = set()
            for h in hit['hsps']:
                covered.update(range(h['query_from'], h['query_to'] + 1))
            in_gene = bool({base(a) for a in hit_accs} & orth)
            rows.append(dict(gene=gene, species=sci, top_hit=','.join(sorted(hit_accs)),
                             hit_title=hit['description'][0]['title'], in_ortholog_gene=in_gene,
                             is_selected_isoform=base(seqs[gene][tax][0]) in {base(a) for a in hit_accs},
                             blast_pct_identity=round(100 * hsp['identity'] / hsp['align_len'], 1),
                             query_coverage_pct=round(100 * len(covered) / search['query_len'], 1),
                             evalue=hsp['evalue'], bit_score=hsp['bit_score']))
            if not in_gene:
                log(f'{gene}: WARNING top {sci} BLAST hit {hit_accs} is not a protein of the NCBI ortholog gene '
                    f'{candidates[gene][tax]["symbol"]} (possible paralog); selected sequence kept, see BLAST table.')
        log(f'{gene}: BLASTp ({report["version"]}, db {report["search_target"]["db"]}) top hit per species checked '
            f'against the ortholog gene (paralog check). Raw output {path.name} accessed {accessed(path)}.')
        setting('5 BLAST QC', 'NCBI BLAST URL API, blastp', report['version'],
                f'db refseq_protein; Entrez query "{entrez}"; hitlist 100; matrix {report["params"].get("matrix")}; '
                f'gap open/extend {report["params"].get("gap_open")}/{report["params"].get("gap_extend")}; '
                f'expect {report["params"].get("expect")}')
    df = pd.DataFrame(rows)
    df.to_csv(D['03_blast'] / 'blast_qc_summary.csv', index=False)
    return df


# ---------- Step 6: MSA ----------
def align(refs, seqs):
    alns = {}
    for gene in GENES:
        fasta = ''.join(raw for tax in SPECIES if tax in seqs[gene]
                        for *_, raw in records((ORIG / f'{gene}_{SPECIES[tax][0].replace(" ", "_")}_{seqs[gene][tax][0]}.fasta').read_text()))
        (D['04_alignments'] / f'{gene}_input.fasta').write_text(fasta, encoding='utf-8')
        path = D['04_alignments'] / f'{gene}_clustalo.aln'

        def run():
            job = get(f'{EBI}/run', data=dict(email=EMAIL, title=f'{gene}_protein', stype='protein',
                                               outfmt='clustal_num', order='input', sequence=fasta))
            while (status := get(f'{EBI}/status/{job}')) in ('QUEUED', 'RUNNING', 'PENDING'):
                time.sleep(5)
            assert status == 'FINISHED', f'Clustal Omega job {job}: {status}'
            log(f'{gene}: Clustal Omega job {job} finished.')
            return get(f'{EBI}/result/{job}/aln-clustal_num')

        text = cached(path, run)
        aln = AlignIO.read(io.StringIO(text), 'clustal')
        AlignIO.write(aln, D['04_alignments'] / f'{gene}_aligned.fasta', 'fasta')
        assert len({len(r) for r in aln}) == 1
        assert str(aln[0].seq).replace('-', '') == refs[gene]['seq'], f'{gene}: human row != reference'
        assert [r.id for r in aln] == [seqs[gene][t][0] for t in SPECIES if t in seqs[gene]]
        longest = max(len(m) for r in aln for m in re.findall('-+', str(r.seq)) or [''])
        alns[gene] = aln
        version = text.splitlines()[0].split(' multiple')[0]
        log(f'{gene}: alignment of {len(aln)} proteins, {aln.get_alignment_length()} columns, longest gap run '
            f'{longest} ({version}, accessed {accessed(path)}).')
    setting('6 Multiple alignment', 'EMBL-EBI Job Dispatcher Clustal Omega REST', version,
            'stype=protein, outfmt=clustal_num, order=input, all other parameters default')
    return alns


# ---------- Step 7: identity and conservation ----------
def conservation(alns, seqs):
    rows, cons = [], {}
    for gene, aln in alns.items():
        tax_of = {seqs[gene][t][0]: t for t in seqs[gene]}
        names = [SPECIES[tax_of[r.id]][1] for r in aln]
        rs = [str(r.seq) for r in aln]
        human, hlen = rs[0], len(rs[0].replace('-', ''))
        for r, name in zip(aln[1:], names[1:]):
            s = str(r.seq)
            same, shared = ident(human, s)
            rows.append(dict(gene=gene, species=SPECIES[tax_of[r.id]][0], accession=r.id, human_length=hlen,
                             species_length=len(s.replace('-', '')), aligned_columns=aln.get_alignment_length(),
                             shared_residue_columns=shared, identical=same, pct_identity=round(100 * same / shared, 2),
                             pct_coverage_of_human=round(100 * shared / hlen, 2),
                             gap_columns=sum((x == '-') != (y == '-') for x, y in zip(human, s)),
                             method='Clustal Omega MSA; identical / columns where both have a residue'))
            assert 0 <= rows[-1]['pct_identity'] <= 100 and 0 <= rows[-1]['pct_coverage_of_human'] <= 100
        matrix = pd.DataFrame([[round(100 * pid(a, b), 2) for b in rs] for a in rs],
                              index=names, columns=names)
        matrix.to_csv(D['05_conservation'] / f'{gene}_pairwise_identity_matrix.csv')
        pos, per = 0, []
        consensus = aln.column_annotations.get('clustal_consensus', ' ' * len(human))
        for col in range(len(human)):
            column = [r[col] for r in rs]
            if column[0] == '-':
                continue
            pos += 1
            per.append(dict(human_position=pos, alignment_column=col + 1, human_aa=column[0],
                            **{n: c for n, c in zip(names, column)},
                            n_other_species_identical=sum(c == column[0] for c in column[1:]),
                            n_other_species=len(column) - 1,
                            fully_conserved=all(c == column[0] for c in column),
                            clustal_symbol=consensus[col].strip() or '(none)'))
        cons[gene] = pd.DataFrame(per)
        cons[gene].to_csv(D['05_conservation'] / f'{gene}_per_position_conservation.csv', index=False)
        log(f'{gene}: {cons[gene].fully_conserved.mean() * 100:.1f}% of human residues identical in all '
            f'{len(rs)} species.')
    df = pd.DataFrame(rows)
    df.to_csv(D['05_conservation'] / 'identity_table.csv', index=False)
    setting('7 Identity/conservation', 'pipeline.py (Python)', 'n/a',
            'identity = identical / columns where both sequences have a residue; coverage = those columns / human length')
    return df, cons


# ---------- Step 8: trees ----------
def trees(alns, seqs):
    rng = random.Random(SEED)
    out = {}
    for gene, aln in alns.items():
        tax_of = {seqs[gene][t][0]: t for t in seqs[gene]}
        names = [f'{SPECIES[tax_of[r.id]][1]} {r.id}' for r in aln]
        rs = [str(r.seq) for r in aln]

        def nj(cols):
            sub = [''.join(s[c] for c in cols) for s in rs]
            lower = [[1 - pid(sub[i], sub[j]) for j in range(i)] + [0]
                     for i in range(len(sub))]
            return DistanceTreeConstructor().nj(DistanceMatrix(names, lower))

        def splits(tree):
            every = frozenset(names)
            found = set()
            for c in tree.find_clades():
                side = frozenset(t.name for t in c.get_terminals())
                side = side if names[0] not in side else every - side
                if 1 < len(side) < len(names) - 1:
                    found.add(side)
            return found

        ncol = aln.get_alignment_length()
        tree = nj(range(ncol))
        counts = Counter()
        for _ in range(BOOTSTRAPS):
            counts.update(splits(nj([rng.randrange(ncol) for _ in range(ncol)])))
        tree.root_at_midpoint()
        every = frozenset(names)
        for c in tree.get_nonterminals():
            side = frozenset(t.name for t in c.get_terminals())
            side = side if names[0] not in side else every - side
            c.confidence = round(100 * counts[side] / BOOTSTRAPS) if 1 < len(side) < len(names) - 1 else None
            c.name = None
        Phylo.write(tree, D['06_phylogeny'] / f'{gene}_nj_bootstrap.nwk', 'newick')
        out[gene] = tree
        log(f'{gene}: NJ tree built; bootstrap support for internal splits: '
            + '; '.join(f'{{{", ".join(sorted(n.split()[0] for n in s))}}} {100 * k / BOOTSTRAPS:.0f}%'
                        for s, k in counts.most_common()) + '.')
    setting('8 Phylogeny', 'Biopython DistanceTreeConstructor (neighbour-joining)', 'Biopython 1.88',
            f'p-distance (gap columns excluded pairwise); {BOOTSTRAPS} bootstrap replicates, seed {SEED}; '
            'midpoint-rooted for display; no substitution-model selection')
    return out


# ---------- Steps 9-10: ClinVar ----------
PLP = {'Pathogenic', 'Likely pathogenic', 'Pathogenic/Likely pathogenic'}
REVIEWED = {'criteria provided, single submitter', 'criteria provided, multiple submitters, no conflicts',
            'reviewed by expert panel', 'practice guideline'}
TITLE = re.compile(r'^(N[MR]_\d+\.\d+)\(([^)]+)\):(c\.\S+)(?: \(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\))?')


def bucket(desc):
    first = desc.split(';')[0].strip()
    if first in PLP:
        return 'Pathogenic/Likely pathogenic'
    if first in {'Benign', 'Likely benign', 'Benign/Likely benign'}:
        return 'Benign/Likely benign'
    if first == 'Uncertain significance':
        return 'Uncertain significance'
    if first.startswith('Conflicting'):
        return 'Conflicting classifications'
    return 'Other'


def clinvar(refs):
    allrows = []
    for gene in GENES:
        term = f'{gene}[gene] AND missense_variant[molecular consequence]'
        es = cached(API / f'clinvar_{gene}_esearch.json', lambda: get(EUTILS + 'esearch.fcgi', params=dict(
            db='clinvar', term=term, retmax=10000, retmode='json', **NCBI)))
        ids = json.loads(es)['esearchresult']['idlist']
        for i in range(0, len(ids), 250):
            batch = json.loads(cached(API / f'clinvar_{gene}_esummary_{i // 250:02d}.json', lambda: get(
                EUTILS + 'esummary.fcgi', data=dict(db='clinvar', id=','.join(ids[i:i + 250]), retmode='json', **NCBI))))['result']
            for uid in batch['uids']:
                v = batch[uid]
                gc = v.get('germline_classification', {})
                m = TITLE.match(v['title'])
                traits = [t['trait_name'] for t in gc.get('trait_set', [])]
                allrows.append(dict(gene=gene, variation_id=int(uid), accession=v.get('accession'), title=v['title'],
                                    transcript=m.group(1) if m else '', title_gene=m.group(2) if m else '',
                                    hgvs_c=m.group(3) if m else '',
                                    hgvs_p=f'p.{m.group(4)}{m.group(5)}{m.group(6)}' if m and m.group(4) else '',
                                    classification=gc.get('description', ''), category=bucket(gc.get('description', '')),
                                    review_status=gc.get('review_status', ''), last_evaluated=gc.get('last_evaluated', ''),
                                    conditions='; '.join(traits),
                                    consequences='; '.join(v.get('molecular_consequence_list', [])),
                                    record_url=f'https://www.ncbi.nlm.nih.gov/clinvar/variation/{uid}/',
                                    date_accessed=accessed(API / f'clinvar_{gene}_esummary_{i // 250:02d}.json')))
        log(f'{gene}: ClinVar query "{term}" returned {len(ids)} records (accessed {accessed(API / f"clinvar_{gene}_esearch.json")}).')
    df = pd.DataFrame(allrows)
    reasons = []
    for r in df.itertuples():
        if r.title_gene != r.gene or not r.hgvs_p or r.hgvs_p.endswith('Ter') or 'missense variant' not in r.consequences:
            reasons.append('not a single-gene missense record for this gene')
        elif r.category != 'Pathogenic/Likely pathogenic':
            reasons.append(f'classification: {r.category}')
        elif r.review_status not in REVIEWED:
            reasons.append('review status below criteria provided')
        elif 'alzheimer' not in r.conditions.lower():
            reasons.append('no Alzheimer condition listed')
        else:
            reasons.append('selected')
    df['selection'] = reasons
    df.to_csv(D['07_variants'] / 'clinvar_all_missense.csv', index=False)
    for gene, g in df.groupby('gene'):
        log(f'{gene}: ClinVar selection outcome: ' + '; '.join(f'{k} {v}' for k, v in g.selection.value_counts().items()) + '.')
    sel = df[df.selection == 'selected'].copy()
    sel.to_csv(D['07_variants'] / 'clinvar_selected.csv', index=False)
    setting('9 ClinVar', 'NCBI E-utilities esearch/esummary, db=clinvar', EUTILS,
            'term "<GENE>[gene] AND missense_variant[molecular consequence]"; selection rule in methods_log.md')
    return df, sel


def map_variants(sel, refs, alns, cons):
    rows = []
    for v in sel.itertuples():
        ref, c = refs[v.gene], cons[v.gene]
        m = TITLE.match(v.title)
        pos, ref_aa, alt_aa = int(m.group(5)), seq1(m.group(4)), seq1(m.group(6))
        base = dict(gene=v.gene, variation_id=v.variation_id, hgvs_p=v.hgvs_p, transcript=v.transcript,
                    classification=v.classification, review_status=v.review_status, residue_position=pos,
                    reference_aa=ref_aa, alternate_aa=alt_aa, record_url=v.record_url)
        same_base = v.transcript.split('.')[0] == ref['nm'].split('.')[0]
        if not same_base:
            rows.append(dict(base, mapping='unmapped', mapping_note=f'transcript {v.transcript} is not reference {ref["nm"]}'))
            continue
        if pos > len(ref['seq']) or ref['seq'][pos - 1] != ref_aa:
            rows.append(dict(base, mapping='unmapped', mapping_note=f'reference residue mismatch at {pos}'))
            continue
        p = c[c.human_position == pos].iloc[0]
        assert p.human_aa == ref_aa
        names = [k for k in SPECIES_NAMES if k in c.columns]
        rows.append(dict(base, mapping='mapped', alignment_column=int(p.alignment_column),
                         residues_in_species=' '.join(f'{n}:{p[n]}' for n in names),
                         n_other_species_identical=int(p.n_other_species_identical),
                         n_other_species=int(p.n_other_species), fully_conserved=bool(p.fully_conserved),
                         clustal_symbol=p.clustal_symbol,
                         mapping_note='' if v.transcript == ref['nm'] else f'version differs from {ref["nm"]}; reference residue verified'))
    df = pd.DataFrame(rows)
    df.to_csv(D['07_variants'] / 'variant_mapping.csv', index=False)
    for gene, g in df.groupby('gene'):
        mapped = g[g.mapping == 'mapped']
        log(f'{gene}: {len(mapped)}/{len(g)} selected variants mapped ({mapped.residue_position.nunique()} distinct '
            f'residues); {int(mapped.fully_conserved.sum())} at residues identical in all species.'
            + (f' Unmapped: {"; ".join(g[g.mapping != "mapped"].mapping_note)}.' if (g.mapping != 'mapped').any() else ''))
    return df


SPECIES_NAMES = [c for _, c in SPECIES.values()]


# ---------- Step 11: UniProt domains and AlphaFold ----------
def structures(refs):
    out = {}
    for gene, ref in refs.items():
        q = f'gene_exact:{gene} AND organism_id:9606 AND reviewed:true'
        hits = json.loads(cached(API / f'uniprot_{gene}_search.json', lambda: get(
            'https://rest.uniprot.org/uniprotkb/search', params=dict(query=q, fields='accession', format='json'))))['results']
        assert len(hits) == 1, f'{gene}: expected one reviewed UniProt entry'
        acc = hits[0]['primaryAccession']
        entry = json.loads(cached(API / f'uniprot_{acc}.json', lambda: get(f'https://rest.uniprot.org/uniprotkb/{acc}.json')))
        feats = pd.DataFrame([dict(type=f['type'], description=f.get('description', ''),
                                   start=f['location']['start']['value'], end=f['location']['end']['value'])
                              for f in entry['features']])
        feats.to_csv(D['08_structures'] / f'{gene}_{acc}_uniprot_features.csv', index=False)
        same = entry['sequence']['value'] == ref['seq']
        log(f'{gene}: UniProt {acc} canonical sequence ' + ('identical to' if same else 'DIFFERS from')
            + f' {ref["np"]}; ' + ('feature coordinates transferred.' if same else 'domain annotations not transferred.'))
        af = json.loads(cached(D['08_structures'] / f'alphafold_{acc}_prediction.json',
                               lambda: get(f'https://alphafold.ebi.ac.uk/api/prediction/{acc}')))[0]
        pdb = cached(D['08_structures'] / Path(af['pdbUrl']).name, lambda: get(af['pdbUrl']))
        plddt = [float(l[60:66]) for l in pdb.splitlines() if l.startswith('ATOM') and l[12:16].strip() == 'CA']
        af_same = af.get('uniprotSequence', af.get('sequence')) == ref['seq']
        log(f'{gene}: AlphaFold DB model {af.get("entryId")} (predicted, not experimental), {len(plddt)} residues, '
            f'mean pLDDT {sum(plddt) / len(plddt):.1f}; sequence ' + ('matches' if af_same else 'DOES NOT match')
            + ' the reference protein. PDB experimental structures not analysed.')
        pd.DataFrame(dict(human_position=range(1, len(plddt) + 1), plddt=plddt)).to_csv(
            D['08_structures'] / f'{gene}_alphafold_plddt.csv', index=False)
        out[gene] = dict(uniprot=acc, features=feats if same else feats.iloc[0:0], plddt=plddt if af_same else [],
                         af_id=af.get('entryId'), af_version=af.get('latestVersion'))
    setting('11 Domains/structure', 'UniProt REST; AlphaFold DB API', 'rest.uniprot.org; alphafold.ebi.ac.uk/api',
            'reviewed human entry; features transferred only if sequence identical; pLDDT from CA B-factor')
    return out


# ---------- Step 12: figures and tables ----------
TRACK_TYPES = {'Transmembrane': '#4c72b0', 'Intramembrane': '#55a868', 'Domain': '#c44e52', 'Active site': '#000000'}


def fig_workflow():
    steps = ['1 Project folders\n& decision log', '2 Human MANE Select\nproteins (NCBI)', '3 Orthologs: 4 mammals\n(NCBI Datasets)',
             '4 FASTA validation\n& metadata sheet', '5 BLASTp QC\n(paralog check)', '6 Clustal Omega\nalignment (EBI)',
             '7 Identity &\nconservation', '8 NJ trees\n(1000 bootstraps)', '9 ClinVar P/LP\nmissense variants',
             '10 Map variants to\nalignment columns', '11 UniProt domains\n& AlphaFold pLDDT', '12 Figures, tables\n& report']
    fig, ax = plt.subplots(figsize=(12, 5.2))
    ax.set_xlim(0, 4)
    ax.set_ylim(0, 3)
    ax.axis('off')
    xy = []
    for i, s in enumerate(steps):
        row, col = divmod(i, 4)
        col = col if row % 2 == 0 else 3 - col  # snake layout
        x, y = col + 0.5, 2.5 - row
        xy.append((x, y))
        ax.text(x, y, s, ha='center', va='center', fontsize=10,
                bbox=dict(boxstyle='round,pad=0.5', fc='#eaf1fb', ec='#2b5797'))
    for (x1, y1), (x2, y2) in zip(xy, xy[1:]):
        dx, dy = x2 - x1, y2 - y1
        ax.annotate('', xy=(x2 - 0.33 * dx, y2 - 0.3 * dy), xytext=(x1 + 0.33 * dx, y1 + 0.3 * dy),
                    arrowprops=dict(arrowstyle='->', color='#2b5797', lw=1.5))
    ax.set_title('Figure 1. Analysis workflow (each step writes to its numbered project folder)', fontsize=12)
    fig.savefig(D['09_figures'] / 'Fig1_workflow.png', dpi=300, bbox_inches='tight')
    plt.close(fig)


def fig_gene(gene, c, varmap, st, num):
    L = len(c)
    fig, axes = plt.subplots(4, 1, figsize=(12, 7.5), sharex=True, gridspec_kw=dict(height_ratios=[3, 1, 2, 2]))
    frac = c.n_other_species_identical / c.n_other_species
    axes[0].bar(c.human_position, frac, width=1, color='#c9d6ea', label='per residue')
    axes[0].plot(c.human_position, frac.rolling(15, center=True, min_periods=1).mean(), color='#2b5797', lw=1.5,
                 label='15-residue moving average')
    axes[0].set_ylabel('Fraction of other\nspecies identical')
    axes[0].set_ylim(0, 1.05)
    axes[0].legend(loc='lower left', fontsize=8, ncol=2)
    f = st['features']
    shown = set()
    for r in f[f.type.isin(TRACK_TYPES) | f.description.str.contains('Amyloid-beta protein 42')].itertuples():
        color = TRACK_TYPES.get(r.type, '#8172b2')
        label = r.type if r.type in TRACK_TYPES else 'Amyloid-beta 42 region'
        if r.type == 'Active site':
            axes[1].plot(r.start, 0.5, 'v', color=color, ms=8, label=None if label in shown else label)
        else:
            axes[1].axvspan(r.start - 0.5, r.end + 0.5, ymin=0.2, ymax=0.8, color=color, alpha=0.7,
                            label=None if label in shown else label)
            if r.type == 'Domain':
                axes[1].text((r.start + r.end) / 2, 0.92, r.description, ha='center', va='bottom', fontsize=7)
        shown.add(label)
    axes[1].set_yticks([])
    axes[1].set_ylim(0, 1.3)
    axes[1].set_ylabel('UniProt\nfeatures')
    axes[1].legend(loc='upper left', bbox_to_anchor=(1.0, 1.0), fontsize=7)
    m = varmap[(varmap.gene == gene) & (varmap.mapping == 'mapped')]
    counts = m.groupby('residue_position').agg(n=('variation_id', 'size'), cons=('fully_conserved', 'first'))
    for p, r in counts.iterrows():
        col = '#b2182b' if r.cons else '#f4a582'
        axes[2].vlines(p, 0, r.n, color=col, lw=1)
        axes[2].plot(p, r.n, 'o', color=col, ms=4)
    axes[2].plot([], [], 'o', color='#b2182b', label='residue identical in all 5 species')
    axes[2].plot([], [], 'o', color='#f4a582', label='residue differs in ≥1 species')
    axes[2].set_ylabel('ClinVar P/LP\nvariants')
    axes[2].set_ylim(0, max(counts.n.max() if len(counts) else 1, 1) + 1)
    axes[2].legend(loc='upper left', fontsize=8)
    if st['plddt']:
        axes[3].plot(range(1, len(st['plddt']) + 1), st['plddt'], color='#444', lw=0.8)
        for y, col in [(90, '#0053d6'), (70, '#65cbf3'), (50, '#ffdb13')]:
            axes[3].axhline(y, color=col, lw=0.8, ls='--')
    axes[3].set_ylabel('AlphaFold pLDDT\n(predicted)')
    axes[3].set_ylim(0, 100)
    axes[3].set_xlabel(f'Human {gene} residue position (reference protein, {L} aa)')
    axes[3].set_xlim(0, L + 1)
    fig.suptitle(f'Figure {num}. {gene}: conservation across human, mouse, rat, dog and cattle, UniProt features, '
                 f'ClinVar P/LP missense variants and AlphaFold confidence', fontsize=10)
    fig.tight_layout()
    fig.savefig(D['09_figures'] / f'Fig{num}_{gene}_conservation_variants.png', dpi=300, bbox_inches='tight')
    plt.close(fig)


def fig_alignment_window(gene, aln, seqs, c, varmap, num, width=60):
    m = varmap[(varmap.gene == gene) & (varmap.mapping == 'mapped')]
    var_cols = set(m.alignment_column - 1)
    n = aln.get_alignment_length()
    start = max(range(max(n - width, 0) + 1), key=lambda s: (sum(s <= v < s + width for v in var_cols), -s))
    tax_of = {seqs[gene][t][0]: t for t in seqs[gene]}
    rows = [(f'{SPECIES[tax_of[r.id]][1]} ({r.id})', str(r.seq)[start:start + width]) for r in aln]
    human = rows[0][1]
    fig, ax = plt.subplots(figsize=(14, 0.45 * len(rows) + 1.6))
    for i, (_, s) in enumerate(rows):
        for j, ch in enumerate(s):
            same = ch == human[j] and ch != '-'
            ax.add_patch(plt.Rectangle((j, -i - 0.5), 1, 1, color='#cfe3cf' if same else '#f6c9c9' if ch != '-' else '#eeeeee'))
            ax.text(j + 0.5, -i, ch, ha='center', va='center', fontsize=8, family='monospace')
    for v in var_cols:
        if start <= v < start + width:
            ax.plot(v - start + 0.5, 0.85, 'v', color='#b2182b', ms=6)
    hpos = c.set_index('alignment_column').human_position
    ticks = [j for j in range(width) if (start + j + 1) in hpos.index and hpos[start + j + 1] % 10 == 0]
    ax.set_xticks([j + 0.5 for j in ticks], [str(hpos[start + j + 1]) for j in ticks])
    ax.set_yticks([-i for i in range(len(rows))], [name for name, _ in rows], fontsize=9)
    ax.set_xlim(0, width)
    ax.set_ylim(-len(rows) + 0.5, 1.3)
    ax.set_xlabel('Human residue position')
    for s in ax.spines.values():
        s.set_visible(False)
    inside = c[(c.alignment_column > start) & (c.alignment_column <= start + width)].human_position
    first, last = inside.min(), inside.max()
    ax.set_title(f'Figure {num}. {gene} Clustal Omega alignment, human residues {first}-{last} (window with most '
                 f'mapped ClinVar P/LP variants, ▼). Green = identical to human, red = different, grey = gap.', fontsize=9)
    fig.savefig(D['09_figures'] / f'Fig{num}_{gene}_alignment_window.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    return first, last


def fig_trees(tree_by_gene, num):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for ax, (gene, tree) in zip(axes, tree_by_gene.items()):
        Phylo.draw(tree, axes=ax, do_show=False,
                   branch_labels=lambda cl: f'{cl.confidence:.0f}' if cl.confidence is not None else None)
        ax.set_title(gene)
        ax.set_xlabel('p-distance (substitutions per site)')
        ax.set_ylabel('')
        ax.set_yticks([])
    fig.suptitle(f'Figure {num}. Neighbour-joining protein trees per gene (p-distance, {BOOTSTRAPS} bootstrap '
                 'replicates; support % on branches; midpoint-rooted for display)', fontsize=10)
    fig.tight_layout()
    fig.savefig(D['09_figures'] / f'Fig{num}_phylogenetic_trees.png', dpi=300, bbox_inches='tight')
    plt.close(fig)


def tables(meta, blast, identity, allvars, varmap):
    settings = pd.DataFrame(SETTINGS)
    settings.to_csv(D['01_metadata'] / 'analysis_settings.csv', index=False)
    with pd.ExcelWriter(D['09_figures'] / 'tables.xlsx') as xw:
        meta.to_excel(xw, sheet_name='A_sequence_metadata', index=False)
        blast.to_excel(xw, sheet_name='BLAST_QC', index=False)
        identity.to_excel(xw, sheet_name='B_pairwise_identity', index=False)
        allvars[allvars.selection == 'selected'].to_excel(xw, sheet_name='ClinVar_selected', index=False)
        varmap.to_excel(xw, sheet_name='C_variant_conservation', index=False)
        allvars.groupby(['gene', 'category']).size().unstack(fill_value=0).to_excel(xw, sheet_name='ClinVar_category_counts')
        settings.to_excel(xw, sheet_name='Analysis_settings', index=False)


def main():
    setup()
    meta = []
    refs = human_refs(meta)
    seqs, candidates = orthologs(refs, meta)
    meta_df = validate(seqs, meta)
    blast = blast_qc(refs, seqs, candidates)
    alns = align(refs, seqs)
    identity, cons = conservation(alns, seqs)
    tree_by_gene = trees(alns, seqs)
    allvars, sel = clinvar(refs)
    varmap = map_variants(sel, refs, alns, cons)
    st = structures(refs)
    fig_workflow()
    windows = {}
    for i, gene in enumerate(GENES):
        fig_gene(gene, cons[gene], varmap, st[gene], f'2{"abc"[i]}')
        windows[gene] = fig_alignment_window(gene, alns[gene], seqs, cons[gene], varmap, f'3{"abc"[i]}')
    fig_trees(tree_by_gene, 4)
    (D['09_figures'] / 'alignment_windows.json').write_text(json.dumps({g: list(map(int, w)) for g, w in windows.items()}))
    setting('12 Figures', 'matplotlib', matplotlib.__version__, '300 dpi PNG')
    tables(meta_df, blast, identity, allvars, varmap)
    log('Figures written to 09_figures/, tables to 09_figures/tables.xlsx, settings to 01_metadata/analysis_settings.csv.')


if __name__ == '__main__':
    main()
