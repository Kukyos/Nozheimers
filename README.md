# Comparative sequence analysis of APP, PSEN1 and PSEN2

This project follows `Alzheimers_Gene_Sequence_Comparative_Analysis_Project_Guide.docx`, the controlling protocol. It compares three early-onset Alzheimer's genes (APP, PSEN1, PSEN2) across human, mouse, rat, dog and cattle. For each gene it measures how conserved the protein is and checks whether curated ClinVar Pathogenic/Likely pathogenic missense variants fall on conserved residues. It is a dry-lab analysis of sequences and database records only. It does not diagnose anything and does not predict anyone's risk.

## Reproduce

Requires Python 3.10+ and an internet connection.

```
pip install -r requirements.txt
python pipeline.py      # Steps 1-12: data, analyses, figures, tables (~5-15 min, mostly BLAST queueing)
python make_report.py   # Steps 13-14: 10_report/report.docx, numbers read from the CSVs
```

- **Caching.** Every download is cached, and its first access date is recorded in `01_metadata/access_dates.json`. A re-run uses the saved files and never overwrites an original.
- **Re-fetching.** Delete a cached file to fetch it again. Expect database records (ClinVar especially) to change over time.
- **Contact email.** `EMAIL` at the top of `pipeline.py` is sent to EBI (required) and NCBI (requested). Change it to your own.

## Folder contents

| Folder | Contents |
|---|---|
| `01_metadata/` | `methods_log.md`: the methods and decision log, with predefined rules first and then everything that was done, including exclusions. `metadata.csv`/`.xlsx`: every sequence considered, with accession.version, URL, date and the include/exclude decision. `analysis_settings.csv`: tools, versions and parameters. `access_dates.json`: first access date of each download. `api_raw/`: raw NCBI, UniProt and ClinVar API responses. |
| `02_raw_sequences/` | `original/`: FASTA exactly as downloaded, one file per selected sequence, plus the human CDS. `candidates/`: every RefSeq isoform downloaded for each ortholog before selection. No cleaned copies were needed. |
| `03_blast/` | Raw BLASTp output (JSON) per gene and `blast_qc_summary.csv` (top hit per species and paralog check). |
| `04_alignments/` | Clustal Omega input FASTA, `.aln` output (clustal_num) and the aligned FASTA, per gene. |
| `05_conservation/` | `identity_table.csv` (Table B), pairwise identity matrices and per-position conservation, per gene. |
| `06_phylogeny/` | Newick trees with bootstrap support, per gene. |
| `07_variants/` | `clinvar_all_missense.csv` (every ClinVar missense record and why it was or wasn't selected), `clinvar_selected.csv` and `variant_mapping.csv` (Table C). |
| `08_structures/` | UniProt feature tables, AlphaFold DB model files and per-residue pLDDT. |
| `09_figures/` | Figures 1-4 (300 dpi PNG) and `tables.xlsx` (all tables, one sheet each). |
| `10_report/` | `report.docx` and the Crossref metadata for every reference. |

## Doing the same steps by hand with free web tools

1. **Human reference.** Open https://www.ncbi.nlm.nih.gov/gene/351 (APP), `/5663` (PSEN1) or `/5664` (PSEN2). Under *NCBI Reference Sequences (RefSeq)*, find the transcript tagged **MANE Select** and its NP_ protein. Open the protein, switch to *FASTA*, then *Send to > File*.
2. **Orthologs.** On the Gene page open *Orthologs*, filter to mouse, rat, dog and cattle, and open each gene. From its RefSeq proteins, choose the isoform that best matches the human protein (`metadata.csv` lists every candidate and its score), then download the FASTA.
3. **BLAST check.** At https://blast.ncbi.nlm.nih.gov/Blast.cgi run *blastp*. Paste the human NP_ accession, set database *RefSeq protein* and organism *Mus musculus, Rattus norvegicus, Canis lupus familiaris, Bos taurus*, and keep the defaults. Check that the top hit for each species is the ortholog, not a paralog (APLP1/APLP2 for APP; PSEN2 for PSEN1 and the reverse).
4. **Alignment.** At https://www.ebi.ac.uk/jdispatcher/msa/clustalo paste the five FASTA records for one gene (human first). Choose output *ClustalW with character counts* and order *input*, keep the rest default, and submit. Repeat for each gene.
5. **Identity.** For each pair, count identical residues and divide by the columns where both sequences have a residue (gap columns excluded). This definition is used everywhere. Clustal Omega's *Percent Identity Matrix* tab uses its own definition, so its numbers can differ slightly.
6. **Trees.** Load the aligned FASTA (`04_alignments/*_aligned.fasta`) into MEGA (https://www.megasoftware.net/). Choose *Phylogeny > Construct Neighbor-Joining Tree* with model *p-distance*, *pairwise deletion* and *1000 bootstrap replicates*. MEGA's *Models > Find Best Protein Models* step was not done here and is a possible extension.
7. **ClinVar.** At https://www.ncbi.nlm.nih.gov/clinvar/ search `APP[gene] AND missense_variant[molecular consequence]`, then apply the selection rule in `01_metadata/methods_log.md`. Check each variant's reference residue against the human protein before mapping it to an alignment column.
8. **Domains and structure.** UniProt: https://www.uniprot.org/uniprotkb/P05067 (APP), P49768 (PSEN1), P49810 (PSEN2); see *Features*. AlphaFold DB: https://alphafold.ebi.ac.uk/entry/P05067 and the same pattern for the others. AlphaFold models are predictions, not experimental structures.

## Not included (each would be a separate extension)

- APOE, which the guide marks as optional.
- Experimental PDB structures.
- Substitution-model selection, because the trees use p-distance.
- A comparison of benign against pathogenic variant positions.

## AI use

The scripts, figures and draft report were produced with the help of an AI assistant (Claude). All data come from the live databases listed above and are traceable to accession.version and access date. Before submitting, check your institution's rules on AI use and disclosure, as section 8 of the guide asks.
