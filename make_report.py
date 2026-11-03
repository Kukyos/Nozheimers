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
TITLE = ('Comparative Genomic and Protein Sequence Analysis of APP, PSEN1 and PSEN2: Evolutionary Conservation and '
         "Contextualization of Curated Alzheimer's Disease-Associated Variants Across Selected Mammalian Species")

REFS = {  # key -> DOI; metadata comes from Crossref, not from memory
    'goate': '10.1038/349704a0', 'sherrington': '10.1038/375754a0', 'levylahad': '10.1126/science.7638622',
    'rogaev': '10.1038/376775a0', 'selkoe': '10.15252/emmm.201606210', 'destrooper': '10.1038/34910',
    'wolfe': '10.1038/19077', 'refseq': '10.1093/nar/gkv1189', 'mane': '10.1038/s41586-022-04558-8',
    'datasets': '10.1038/s41597-024-03571-y', 'blast': '10.1186/1471-2105-10-421', 'clustalo': '10.1038/msb.2011.75',
    'jdispatcher': '10.1093/nar/gkae241', 'blosum': '10.1073/pnas.89.22.10915', 'nj': '10.1093/oxfordjournals.molbev.a040454',
    'bootstrap': '10.1111/j.1558-5646.1985.tb00420.x', 'biopython': '10.1093/bioinformatics/btp163',
    'clinvar': '10.1093/nar/gkx1153', 'uniprot': '10.1093/nar/gkae1010', 'alphafold': '10.1038/s41586-021-03819-2',
    'afdb': '10.1093/nar/gkab1061', 'matplotlib': '10.1109/MCSE.2007.55',
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
