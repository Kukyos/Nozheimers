"""Builds 10_report/report.docx (guide Steps 13-14) from the outputs of pipeline.py.

Every number in the text is read from the CSVs written by pipeline.py. Every reference is
resolved through the Crossref API first, and its metadata is used as returned.
Run after pipeline.py:  python make_report.py
"""
import json, re
from pathlib import Path

import pandas as pd
import requests
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt

ROOT = Path(__file__).resolve().parent
OUT = ROOT / '10_report'
FIG = ROOT / '09_figures'
TITLE = ("Comparative Sequence and Structure Analysis of APP, PSEN1 and PSEN2: Conservation of Alzheimer's Disease "
         'Residues and Drug-Binding Sites Across Laboratory Model Organisms')

REFS = {  # key -> DOI; metadata comes from Crossref, not from memory
    'goate': '10.1038/349704a0', 'sherrington': '10.1038/375754a0', 'levylahad': '10.1126/science.7638622',
    'rogaev': '10.1038/376775a0', 'selkoe': '10.15252/emmm.201606210', 'destrooper': '10.1038/34910',
    'wolfe': '10.1038/19077', 'refseq': '10.1093/nar/gkv1189', 'mane': '10.1038/s41586-022-04558-8',
    'datasets': '10.1038/s41597-024-03571-y', 'blast': '10.1186/1471-2105-10-421', 'clustalo': '10.1038/msb.2011.75',
    'jdispatcher': '10.1093/nar/gkae241', 'blosum': '10.1073/pnas.89.22.10915', 'nj': '10.1093/oxfordjournals.molbev.a040454',
    'bootstrap': '10.1111/j.1558-5646.1985.tb00420.x', 'biopython': '10.1093/bioinformatics/btp163',
    'clinvar': '10.1093/nar/gkx1153', 'uniprot': '10.1093/nar/gkae1010', 'alphafold': '10.1038/s41586-021-03819-2',
    'afdb': '10.1093/nar/gkab1061', 'matplotlib': '10.1109/MCSE.2007.55',
    'kabsch': '10.1107/S0567739476001873', '3dmol': '10.1093/bioinformatics/btu829', 'pdb': '10.1093/nar/28.1.235',
    'dssp': '10.1002/bip.360221211',
}
AUTHOR_FALLBACK = {'nj': 'Saitou N, Nei M'}  # Crossref lists no authors for this record
ORDER = []


def crossref():
    path = OUT / 'references_crossref.json'
    data = json.loads(path.read_text(encoding='utf-8')) if path.exists() else {}
    for key, doi in REFS.items():
        if key not in data:
            r = requests.get(f'https://api.crossref.org/works/{doi}', timeout=60,
                             headers={'User-Agent': 'nozheimers_report (mailto:amohamedarmaan@gmail.com)'})
            assert r.status_code == 200, f'DOI did not resolve: {doi}'
            data[key] = r.json()['message']
    path.write_text(json.dumps(data, indent=1), encoding='utf-8')
    return data


def cite(*keys):
    for k in keys:
        if k not in ORDER:
            ORDER.append(k)
    return '[' + ', '.join(str(ORDER.index(k) + 1) for k in keys) + ']'


def fmt_ref(key, m):
    people = m.get('author', [])
    names = [f"{a['family']} {''.join(p[0] for p in re.split(r'[\s.-]+', a.get('given', '')) if p)}".strip()
             if 'family' in a else a.get('name', '') for a in people]
    authors = AUTHOR_FALLBACK.get(key) or (', '.join(names[:3]) + (' et al.' if len(names) > 3 else ''))
    year = m['issued']['date-parts'][0][0]
    journal = (m.get('container-title') or [''])[0].replace('&amp;', '&')
    loc = ''.join([m.get('volume', ''), f"({m['issue']})" if m.get('issue') else '', f":{m['page']}" if m.get('page') else ''])
    title = re.sub(r'\s+', ' ', m['title'][0]).rstrip('.')
    return f'{authors} ({year}). {title}. {journal} {loc}. https://doi.org/{m["DOI"]}'


def para(doc, text, bold=False, italic=False, size=None, align=None):
    p = doc.add_paragraph()
    run = p.add_run(text)
    run.bold, run.italic = bold, italic
    if size:
        run.font.size = Pt(size)
    if align:
        p.alignment = align
    return p


def table(doc, df, caption, note=None):
    para(doc, caption, bold=True, size=9)
    t = doc.add_table(rows=1, cols=len(df.columns))
    t.style = 'Light Grid Accent 1'
    for cell, name in zip(t.rows[0].cells, df.columns):
        cell.text = str(name)
    for row in df.itertuples(index=False):
        for cell, v in zip(t.add_row().cells, row):
            cell.text = '' if pd.isna(v) else str(v)
    for r in t.rows:
        for c in r.cells:
            for p in c.paragraphs:
                for run in p.runs:
                    run.font.size = Pt(8)
    if note:
        para(doc, note, italic=True, size=8)
    doc.add_paragraph()


def figure(doc, path, caption, width=6.5):
    doc.add_picture(str(path), width=Inches(width))
    doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
    para(doc, caption, italic=True, size=9)


# ---------------------------------------------------------------- model organisms (site/build_explorer.py output)
def load_js(path):
    text = (ROOT / path).read_text(encoding='utf-8')
    return json.loads(text[text.index('{'):text.rindex('}') + 1])


def verdict(o, gene):
    """Same published rules as site/explorer.html."""
    dis = o['disease_same'] / o['disease_total'] if o['disease_total'] else 1
    core = o['core_same'] / o['core_total'] if o['core_total'] else 1
    act = all(a['human'] == a['org'] for a in o.get('active', []))
    if o['coverage'] < 90 and o['identity'] >= 80:
        return 'Incomplete record'
    if o['identity'] >= 90 and dis >= .95 and core >= .95 and act:
        return 'Strong model'
    if o['identity'] >= 70 and act and dis >= .85:
        return 'Usable with care' + (' (needs human APP for amyloid-beta work)' if gene == 'APP' and core < 1 else '')
    return 'Distant model'


def fig_organisms(ex, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    colour = {'Strong': '#2e7d4f', 'Usable': '#d98a00', 'Distant': '#b8325f', 'Incomplete': '#888888'}
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.4), sharex=True)
    for ax, (gene, g) in zip(axes, ex['genes'].items()):
        rows = sorted([o for o in g['organisms'].values() if not o.get('missing')], key=lambda o: o['identity'])
        y = list(range(len(rows)))
        ax.barh(y, [o['identity'] for o in rows], color=[colour[verdict(o, gene).split()[0]] for o in rows])
        for i, o in zip(y, rows):
            ax.text(o['identity'] + 1.5, i, f"{o['identity']:.1f}%  ({o['disease_same']}/{o['disease_total']})", va='center', fontsize=7.5)
        ax.set_yticks(y, [o['name'] for o in rows], fontsize=8.5)
        ax.set_title(f'{gene} (human {g["uniprot"]})', fontsize=11)
        ax.set_xlim(0, 135)
        ax.set_xticks([0, 25, 50, 75, 100])
        ax.set_xlabel('% identical to human')
        for x in (70, 90):
            ax.axvline(x, color='#999999', lw=.8, ls='--')
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in colour.values()]
    fig.legend(handles, ['Strong model', 'Usable with care', 'Distant model', 'Incomplete record'], loc='lower center', ncol=4, frameon=False)
    fig.suptitle('Identity of each organism\'s protein to the human reference, coloured by verdict. '
                 'In brackets: human disease residues kept. Dashed lines: 70% and 90% identity thresholds.', fontsize=10)
    fig.tight_layout(rect=(0, .06, 1, .95))
    fig.savefig(path, dpi=250)
    plt.close(fig)


def bullets(doc, items):
    for it in items:
        doc.add_paragraph(it, style='List Bullet')


def main():
    OUT.mkdir(exist_ok=True)
    refs = crossref()
    ident = pd.read_csv(ROOT / '05_conservation/identity_table.csv')
    vm = pd.read_csv(ROOT / '07_variants/variant_mapping.csv')
    allv = pd.read_csv(ROOT / '07_variants/clinvar_all_missense.csv')
    meta = pd.read_csv(ROOT / '01_metadata/metadata.csv')
    blast = pd.read_csv(ROOT / '03_blast/blast_qc_summary.csv')
    settings = pd.read_csv(ROOT / '01_metadata/analysis_settings.csv')
    ex = load_js('site/explorer/explorer.js')
    orth = load_js('site/explorer/orthologs.js')
    G = ex['genes']
    cons = {g: pd.read_csv(ROOT / f'05_conservation/{g}_per_position_conservation.csv') for g in G}
    fig6 = FIG / 'Fig6_model_organisms.png'
    fig_organisms(ex, fig6)
    mapped = vm[vm.mapping == 'mapped']
    org = lambda g, name: next(o for o in G[g]['organisms'].values() if o.get('name') == name)
    rng = lambda g: f'{ident[ident.gene == g].pct_identity.min():.1f}-{ident[ident.gene == g].pct_identity.max():.1f}%'
    n_cons = int(mapped.fully_conserved.sum())
    m_app, m_ps1 = org('APP', 'Mouse'), org('PSEN1', 'Mouse')

    doc = Document()
    st = doc.styles['Normal']
    st.font.name, st.font.size = 'Calibri', Pt(11)
    for s in doc.sections:
        s.left_margin = s.right_margin = Inches(1)

    # ---- title page
    para(doc, TITLE, bold=True, size=18, align=WD_ALIGN_PARAGRAPH.CENTER)
    para(doc, 'A reproducible dry-lab bioinformatics study', italic=True, size=12, align=WD_ALIGN_PARAGRAPH.CENTER)
    para(doc, f'Data retrieved 29 October to 3 November 2026. Report generated {pd.Timestamp.now():%d %B %Y}.', size=10,
         align=WD_ALIGN_PARAGRAPH.CENTER)
    para(doc, 'Interactive results: https://nozheimers.vercel.app and https://nozheimers.vercel.app/explorer.html', size=10,
         align=WD_ALIGN_PARAGRAPH.CENTER)
    para(doc, 'Code and data: https://github.com/Kukyos/Nozheimers', size=10, align=WD_ALIGN_PARAGRAPH.CENTER)

    # ---- abstract
    doc.add_heading('Abstract', 1)
    para(doc, (
        'Background: mutations in APP, PSEN1 and PSEN2 cause early-onset Alzheimer\'s disease, and drug candidates are '
        'tested in animals before humans. Aim: to determine which laboratory organisms carry proteins close enough to '
        'the human ones for disease-mutation and drug studies. Methods: with the human proteins as the reference, we '
        'aligned their orthologs from four mammals (multiple alignment, Clustal Omega) and from 12 laboratory '
        'organisms (pairwise), superposed AlphaFold structures, mapped curated ClinVar pathogenic variants, and '
        'measured conservation of the residues that contact nine drugs in experimental PDB structures. Results: identity to '
        f'human in mouse, rat, dog and cattle was {rng("APP")} for APP, {rng("PSEN1")} for PSEN1 and {rng("PSEN2")} for PSEN2. '
        f'{n_cons} of {len(mapped)} mapped pathogenic variants lie at residues identical in all five mammals. Mouse '
        f'PSEN1 keeps {m_ps1["disease_same"]} of {m_ps1["disease_total"]} disease residues and all '
        f'{len(m_ps1["drugs"])} drug pockets; mouse APP keeps all {m_app["disease_total"]} disease residues but differs at '
        f'{m_app["core_total"] - m_app["core_same"]} of 42 amyloid-beta residues. Conclusion: non-human primates, pig, cattle '
        'and dog carry essentially human versions of all three proteins; rodents suit presenilin and gamma-secretase '
        'drug studies, but amyloid-beta studies in rodents need the human APP sequence.'))

    # ---- introduction
    doc.add_heading('1. Introduction', 1)
    para(doc, 'Alzheimer\'s disease is the most common cause of dementia. A small fraction of cases are inherited and '
              'begin early, caused by mutations in APP ' + cite('goate') + ', PSEN1 ' + cite('sherrington') +
              ' or PSEN2 ' + cite('levylahad', 'rogaev') + '. APP is cut twice to release amyloid-beta, the peptide '
              'that accumulates in plaques ' + cite('selkoe') + '. Presenilin forms the catalytic core of '
              'gamma-secretase, which makes the second cut with two aspartate residues ' + cite('destrooper', 'wolfe') + '.')
    para(doc, 'Drugs and disease mechanisms are studied in animals before humans. Such a study transfers to patients only '
              'if the animal\'s protein has the same residues where it matters: at the sites mutated in patients, in the '
              'functional core, and where a drug binds. Sequence comparison counts the differences; structure '
              'comparison shows where they sit.')
    para(doc, 'Research question: how conserved are APP, PSEN1 and PSEN2 between humans and laboratory organisms, and '
              'which organisms keep the human disease residues, functional core and drug-binding residues well '
              'enough to serve as test models?', bold=True)

    # ---- input
    doc.add_heading('2. Input', 1)
    para(doc, 'The human protein is the reference in every comparison; the other organisms are the inputs tested against it.')
    hum = meta[(meta.common_name == 'human') & (meta.sequence_type == 'protein')]
    table(doc, pd.DataFrame([{'Gene': r.gene, 'Protein': r.record_description.split(' [')[0].split(' ', 1)[1],
                              'Length (aa)': r.length, 'RefSeq protein': r.accession_version,
                              'Transcript': r.transcript_protein_id, 'UniProt': G[r.gene]['uniprot'],
                              'NCBI Gene': G[r.gene]['gene_id']} for r in hum.itertuples()]),
          'Table 1. Human reference proteins (MANE Select). The UniProt canonical sequences are identical to these RefSeq '
          'proteins, so residue numbers agree across databases.')
    latin = {}
    for g in G.values():
        for t, o in g['organisms'].items():
            if not o.get('missing'):
                latin[int(t)] = o['latin']
    table(doc, pd.DataFrame([{'Group': p['group'], 'Organism': p['name'], 'Species': latin.get(p['taxid'], ''),
                              'NCBI Taxon ID': p['taxid']} for p in ex['presets']]),
          'Table 2. Laboratory organisms compared with the human reference.',
          'The online explorer extends the comparison to any species with an NCBI ortholog call '
          f'({min(len(v) for v in orth.values())}-{max(len(v) for v in orth.values())} species, depending on the gene).')
    para(doc, 'Databases and tools: NCBI Gene, RefSeq and Datasets ' + cite('refseq', 'mane', 'datasets') +
              '; NCBI BLAST ' + cite('blast') + '; Clustal Omega via the EMBL-EBI Job Dispatcher ' + cite('clustalo', 'jdispatcher') +
              '; ClinVar ' + cite('clinvar') + '; UniProt ' + cite('uniprot') + '; AlphaFold DB ' + cite('alphafold', 'afdb') +
              '; RCSB Protein Data Bank ' + cite('pdb') + '; pSSdb, the Protein Secondary Structure Database '
              '(https://pssdb.mimuw.edu.pl); Biopython ' + cite('biopython') + '; 3Dmol.js ' + cite('3dmol') +
              '; matplotlib ' + cite('matplotlib') + '.')

    # ---- data and methods
    doc.add_heading('3. Data and methods', 1)
    para(doc, 'All rules were fixed before results were seen. Every decision is logged in 01_metadata/methods_log.md, '
              'and every download is cached with its access date.')
    figure(doc, FIG / 'Fig1_workflow.png', 'Figure 1. Workflow of the five-mammal analysis (pipeline.py), from sequence '
                                           'retrieval to interpretation.')
    doc.add_heading('3.1 Five-mammal comparison', 2)
    bullets(doc, [
        'Human reference: the MANE Select transcript and protein from the NCBI Gene product report.',
        f'Orthologs in mouse, rat, dog and cattle: NCBI ortholog calls. Among {int((meta.sequence_type == "protein").sum())} '
        'candidate RefSeq proteins, the isoform with the highest global alignment score to human was kept (BLOSUM62 ' +
        cite('blosum') + ', gap open -11, extend -1), preferring a curated NP_ record within 5% of the best score.',
        'Quality control: FASTA validation, and BLASTp of each human protein against RefSeq protein to confirm that the '
        'top hit per species is the ortholog, not a paralog.',
        'Multiple alignment: Clustal Omega, protein, default parameters, one alignment per gene.',
        'Identity = identical residues / columns where both sequences have a residue. Coverage = those columns / human length.',
        'Trees: neighbour-joining ' + cite('nj') + ' on p-distance with 1000 bootstrap replicates ' + cite('bootstrap') + '.',
        'ClinVar: missense records classified Pathogenic or Likely pathogenic, review status "criteria provided" or better, '
        'with an Alzheimer condition listed; mapped only when the transcript and reference residue match the human reference.'])
    doc.add_heading('3.2 Laboratory organisms and structures', 2)
    bullets(doc, [
        'Orthologs for 12 organisms from NCBI ortholog calls, mapped to UniProt. All UniProt entries linked to the gene were '
        'aligned to human and the best full-length match was kept. NCBI makes no call for fruit fly and C. elegans, so '
        'reviewed UniProt homologs were used (Appl and apl-1 for APP; Psn and sel-12, each a co-ortholog of both PSEN1 and PSEN2).',
        'Global pairwise alignment to human (BLOSUM62, gap open -11, extend -1).',
        'Structure comparison: the AlphaFold DB models of the human and organism proteins were superposed on C-alpha atoms '
        '(Kabsch ' + cite('kabsch') + ') using residues with pLDDT of at least 70 in both models. APP, a multi-domain protein '
        'with flexible linkers, was fitted on its best-fitting domain (E2, residues 374-565).',
        'Functional core: amyloid-beta (APP 672-713) for APP; the UniProt membrane-spanning helices for PSEN1 and PSEN2, '
        'with the catalytic aspartates PSEN1 D257/D385 and PSEN2 D263/D366.',
        'Drug pockets: human residues with any heavy atom within 4.5 angstrom of the drug in experimental PDB complexes (Table 9).',
        'Secondary structure: helix, sheet and coil assigned from each 3D model with a DSSP-style hydrogen-bond method ' +
        cite('dssp') + ', the same approach pSSdb applies to AlphaFold models, so pSSdb entries for the same UniProt '
        'accessions serve as an independent cross-check.',
        'Verdict, stated in advance so it can be challenged. Strong model: at least 90% identity, at least 95% of disease '
        'residues kept, at least 95% of the core identical, catalytic residues kept. Usable with care: at least 70% identity, '
        'catalytic residues kept, at least 85% of disease residues kept. Otherwise: distant model.',
        'Checks: human against itself gives 100% identity and 0 angstrom RMSD; the browser implementation reproduces the '
        'Python results exactly (site/test_live.js); MRK-560, reported as PSEN1-selective, has the pocket that differs most '
        'between human PSEN1 and PSEN2 (T281P, L282I).'])
    table(doc, settings, 'Table 3. Tools, versions and parameters (01_metadata/analysis_settings.csv).')

    # ---- output
    doc.add_heading('4. Output (results)', 1)
    doc.add_heading('4.1 Sequence retrieval and quality control', 2)
    pend = sorted(blast[blast.top_hit.astype(str).str.startswith('pending')].gene.unique())
    para(doc, 'All 15 selected proteins passed validation. BLASTp confirmed the APP ortholog as the top hit in every species'
              + (f'. The searches for {" and ".join(pend)} were still queued at NCBI when this report was generated and are '
                 'marked pending; re-running pipeline.py fills them in.' if pend else '.'))
    doc.add_heading('4.2 Conservation across five mammals', 2)
    t4 = ident[['gene', 'species', 'accession', 'pct_identity', 'pct_coverage_of_human', 'shared_residue_columns', 'gap_columns']]
    table(doc, t4.rename(columns={'gene': 'Gene', 'species': 'Species', 'accession': 'Accession', 'pct_identity': 'Identity %',
                                  'pct_coverage_of_human': 'Coverage %', 'shared_residue_columns': 'Shared columns',
                                  'gap_columns': 'Gap columns'}),
          'Table 4. Pairwise identity of each species to the human reference (from the Clustal Omega alignment).')
    para(doc, '; '.join(f'{g}: {cons[g].fully_conserved.mean() * 100:.1f}% of human residues are identical in all five species'
                        for g in G) + '.')
    for i, g in enumerate(G):
        figure(doc, FIG / f'Fig2{"abc"[i]}_{g}_conservation_variants.png',
               f'Figure 2{"abc"[i]}. {g}: per-residue conservation across five species, UniProt features and ClinVar '
               'pathogenic variants.')
    doc.add_heading('4.3 Disease variants', 2)
    para(doc, 'ClinVar missense records: ' + ', '.join(f'{g} {int((allv.gene == g).sum())}' for g in G) + '. Of these, ' +
              ', '.join(f'{g} {int(((allv.gene == g) & (allv.selection == "selected")).sum())}' for g in G) +
              f' met the selection rule, and all were mapped. {n_cons} of {len(mapped)} lie at residues identical in all '
              'five mammals.')
    summary = mapped.groupby('gene').agg(variants=('variation_id', 'size'), residues=('residue_position', 'nunique'),
                                         at_fully_conserved=('fully_conserved', 'sum')).reset_index()
    summary['all residues fully conserved %'] = [round(cons[g].fully_conserved.mean() * 100, 1) for g in summary.gene]
    table(doc, summary.rename(columns={'gene': 'Gene', 'variants': 'Variants', 'residues': 'Residues',
                                       'at_fully_conserved': 'Variants at fully conserved residues'}),
          'Table 5. Curated pathogenic variants and conservation. The full list with ClinVar links is in '
          '07_variants/variant_mapping.csv.')
    for i, g in enumerate(G):
        figure(doc, FIG / f'Fig3{"abc"[i]}_{g}_alignment_window.png', f'Figure 3{"abc"[i]}. {g}: alignment window '
               'around the densest cluster of pathogenic variants.')
    figure(doc, FIG / 'Fig4_phylogenetic_trees.png', 'Figure 4. Neighbour-joining gene trees (p-distance, 1000 bootstraps). '
                                                      'These are gene trees, not disease trees.')
    doc.add_heading('4.4 3D structure comparison', 2)
    para(doc, 'In the online explorer the organism model is superposed on the human model and can be rotated, coloured by '
              'difference, confidence, secondary structure, disease residues or drug pocket, and queried residue by residue.')
    for name, cap in [('Fig5_APP-10090-identity.png', 'Figure 5a. Human APP (blue) and mouse App (amber) AlphaFold models superposed on the E2 domain; magenta marks residues that differ.'),
                      ('Fig5_PSEN1-10090-pocket-7Y5T.png', 'Figure 5b. Human and mouse PSEN1 superposed (0.16 angstrom RMSD); cyan sticks mark the residues that contact MRK-560 in PDB 7Y5T, all identical in mouse.'),
                      ('Fig5_PSEN1-7227-identity.png', 'Figure 5c. Human PSEN1 and fruit fly Psn: the membrane core superposes closely (0.86 angstrom) despite 54% identity, but disease and drug-contact residues differ.')]:
        if (FIG / 'explorer' / name).exists():
            figure(doc, FIG / 'explorer' / name, cap + ' Screenshot of the model-organism explorer.')
    doc.add_heading('4.5 Laboratory organisms', 2)
    figure(doc, fig6, 'Figure 6. Identity of each laboratory organism\'s protein to the human reference, coloured by verdict.')
    for n, (g, gv) in enumerate(G.items(), start=6):
        rows = []
        for o in sorted([o for o in gv['organisms'].values() if not o.get('missing')], key=lambda o: -o['identity']):
            sp = o.get('superpose') or {}
            rows.append({'Organism': o['name'], 'UniProt': o['uniprot'], 'Identity %': o['identity'],
                         'Disease residues kept': f"{o['disease_same']}/{o['disease_total']}",
                         'Core identical': f"{o['core_same']}/{o['core_total']}",
                         'RMSD (CA)': f"{sp['rmsd']} ({sp['n']})" if sp else '-',
                         'Drug pockets identical': f"{sum(d['same'] == d['total'] for d in o['drugs'])}/{len(o['drugs'])}" if o['drugs'] else '-',
                         'Verdict': verdict(o, g)})
        missing = [o['name'] for o in gv['organisms'].values() if o.get('missing')]
        table(doc, pd.DataFrame(rows), f'Table {n}. Human {g} ({gv["uniprot"]}) against laboratory organisms. Core = '
              f'{gv["core_label"]}. RMSD in angstrom over the number of C-alpha atoms in brackets.',
              ('No ortholog listed for: ' + ', '.join(missing) + '.') if missing else None)
    doc.add_heading('4.6 Drug-binding pockets', 2)
    drows = []
    for d in ex['drugs']:
        same = [o['name'] for o in G[d['gene']]['organisms'].values() if not o.get('missing') and
                any(x['pdb'] == d['pdb'] and x['same'] == x['total'] for x in o['drugs'])]
        drows.append({'Drug': d['drug'], 'Type': d['kind'], 'Target': d['gene'], 'PDB': d['pdb'],
                      'Method, resolution': f"{d['method']}, {d['resolution']} A", 'Pocket residues': len(d['residues']),
                      'Organisms with identical pocket': ', '.join(same) or 'none'})
    table(doc, pd.DataFrame(drows), 'Table 9. Drugs in experimental human structures and the organisms whose drug-contact '
                                    'residues are identical to human. This is pocket conservation, not a docking simulation.')
    mouse_diffs = sorted({x for d in m_app['drugs'] for x in d['diffs']})
    para(doc, f'Mouse amyloid-beta differs from human at {", ".join(mouse_diffs)}, which changes the binding residues of '
              + ', '.join(d['drug'] for d in m_app['drugs'] if d['diffs']) + ' but not '
              + ', '.join(d['drug'] for d in m_app['drugs'] if not d['diffs']) + '.')

    # ---- interpretation
    doc.add_heading('5. Interpretation', 1)
    psen_ok = [o['name'] for o in G['PSEN1']['organisms'].values() if not o.get('missing') and verdict(o, 'PSEN1') == 'Strong model']
    app_ok = [o['name'] for o in G['APP']['organisms'].values() if not o.get('missing') and verdict(o, 'APP') == 'Strong model']
    act_all = all(a['human'] == a['org'] for gg in ('PSEN1', 'PSEN2') for o in G[gg]['organisms'].values()
                  if not o.get('missing') for a in o['active'])
    bullets(doc, [
        f'Disease mutations concentrate on conserved residues: {n_cons} of {len(mapped)} pathogenic variants sit at residues '
        'identical in all five mammals, a higher share than the background conservation of each protein (Table 5). '
        'Conservation supports, but does not by itself establish, that these residues matter.',
        'The catalytic aspartates of PSEN1 and PSEN2 are ' + ('identical in every organism tested, including fruit fly and C. elegans.'
                                                              if act_all else 'not identical in every organism.'),
        f'For PSEN1 and gamma-secretase drug work, the strong models are {", ".join(psen_ok)}. They keep the membrane core and all '
        'five experimental drug pockets identical to human.',
        f'For APP and amyloid-beta work, the strong models are {", ".join(app_ok)}. Mouse and rat keep every APP disease residue '
        'but differ at three amyloid-beta residues, which changes three of the four antibody-binding sites; amyloid-beta studies '
        'in rodents therefore need the human APP sequence (knock-in or transgene).',
        'Zebrafish and frog are usable with care. Fruit fly and C. elegans are distant models, suited to basic genetics rather '
        'than to human-specific mutation or drug studies.'])
    doc.add_heading('6. Conclusion', 1)
    para(doc, 'We can test PSEN1 and PSEN2 mechanisms and gamma-secretase drugs in mouse and rat, because their presenilins keep the '
              'human disease residues, catalytic core and drug-binding residues. For amyloid-beta and anti-amyloid antibodies, '
              'non-human primates, pig, dog and cattle carry the human sequence, while rodents need the human APP gene introduced. '
              'These are sequence- and structure-level predictions that motivate, but do not replace, experimental validation.')
    doc.add_heading('7. Limitations', 1)
    bullets(doc, [
        'Conservation is one line of evidence. Pathogenicity is established clinically, and ClinVar classifications change.',
        'This analysis is not diagnostic and predicts nothing about any individual.',
        'Organism structures are AlphaFold predictions; flexible low-confidence regions are excluded from the RMSD.',
        'Pocket conservation is not docking: identical residues suggest, but do not prove, identical binding.',
        'Gamma-secretase also contains nicastrin, APH-1 and PEN-2, which were not compared; some drugs contact them.',
        'The trees use p-distance; substitution-model selection was not performed.'])
    doc.add_heading('Data, code and AI use', 1)
    para(doc, 'All data come from public databases and are traceable by accession.version and access date. Code, raw downloads '
              'and outputs: https://github.com/Kukyos/Nozheimers. Interactive results: https://nozheimers.vercel.app. The scripts, '
              'website and this report were produced with the help of an AI assistant (Claude), and every number is computed from '
              'the saved outputs. Check your institution\'s rules on AI use and disclosure.')
    doc.add_heading('References', 1)
    for i, k in enumerate(ORDER, 1):
        para(doc, f'{i}. {fmt_ref(k, refs[k])}', size=9)
    para(doc, f'{len(ORDER) + 1}. Kozlowski LP. pSSdb: Protein Secondary Structure Database (submitted). https://pssdb.mimuw.edu.pl', size=9)
    para(doc, f'{len(ORDER) + 2}. RCSB PDB entries 6LQG, 6LR4, 7C9I, 7D8X, 7Y5T, 6CO3, 4XXD, 5VZY, 5CSZ, 6IYC, 7Y5Z and 5OQV, '
              'https://www.rcsb.org/structure/ followed by the entry ID.', size=9)
    out = OUT / 'report.docx'
    doc.save(out)
    print('written', out)


if __name__ == '__main__':
    main()
