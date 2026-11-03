# Methods and decision log

Run started 2026-11-03 19:39. Rules are fixed here before any retrieval; later entries record what was done.

## Predefined rules
- Genes: APP (Gene ID 351), PSEN1 (5663), PSEN2 (5664), each analysed separately. APOE not included (optional extension in the guide).
- Species: Homo sapiens (human, taxid 9606), Mus musculus (mouse, taxid 10090), Rattus norvegicus (rat, taxid 10116), Canis lupus familiaris (dog, taxid 9615), Bos taurus (cattle, taxid 9913). A species is excluded for a gene if NCBI lists no ortholog; exclusions are logged.
- Human reference: the MANE Select transcript and its protein, read from the NCBI Gene product report.
- Ortholog isoform rule: among all RefSeq proteins of the NCBI-designated ortholog gene, take the one with the highest global alignment score to the human reference (BLOSUM62, gap open -11, extend -1); a curated NP_ record is preferred over a predicted XP_ record when its score is within 5% of the best.
- Sequence type: protein for all alignment, identity, phylogeny and variant analyses. Human CDS is saved for reference only and not mixed into protein analyses.
- Percent identity: identical residues / alignment columns where both sequences have a residue (gap columns excluded). Coverage: those shared columns / human protein length. The same definition gives the p-distance (1 - identity) used for the trees.
- Trees: neighbour-joining on p-distance, 1000 bootstrap replicates (columns resampled, seed 1), midpoint-rooted for display only.
- ClinVar selection rule: missense single-gene records for the gene whose germline classification is Pathogenic, Likely pathogenic or Pathogenic/Likely pathogenic; review status at least "criteria provided, single submitter"; at least one condition name containing "Alzheimer". All other records are kept in the full table but not mapped.
- Variant mapping rule: map only when the ClinVar HGVS transcript is the human reference transcript and the reference amino acid matches the reference protein at that position; otherwise label "unmapped".

## Work log
- APP: 11 RefSeq transcripts listed; selected MANE Select NM_000484.4 -> NP_000475.1 (770 aa). Saved APP_Homo_sapiens_NP_000475.1.fasta and APP_Homo_sapiens_NM_000484.4_CDS.fasta.
- PSEN1: 18 RefSeq transcripts listed; selected MANE Select NM_000021.4 -> NP_000012.1 (467 aa). Saved PSEN1_Homo_sapiens_NP_000012.1.fasta and PSEN1_Homo_sapiens_NM_000021.4_CDS.fasta.
- PSEN2: 18 RefSeq transcripts listed; selected MANE Select NM_000447.3 -> NP_000438.2 (448 aa). Saved PSEN2_Homo_sapiens_NP_000438.2.fasta and PSEN2_Homo_sapiens_NM_000447.3_CDS.fasta.
- APP: Mus musculus ortholog App (Gene ID 11820) has 8 RefSeq proteins; selected NP_001185752.1 (770 aa, score 3923; best candidate score 3923).
- APP: Rattus norvegicus ortholog App (Gene ID 54226) has 11 RefSeq proteins; selected NP_062161.1 (770 aa, score 3937; best candidate score 3937).
- APP: Canis lupus familiaris ortholog APP (Gene ID 403407) has 7 RefSeq proteins; selected NP_001006601.2 (770 aa, score 3947; best candidate score 3947).
- APP: Bos taurus ortholog APP (Gene ID 280722) has 7 RefSeq proteins; selected XP_005201194.1 (770 aa, score 3974; best candidate score 3974); FLAGS: predicted (XP_) record.
- PSEN1: Mus musculus ortholog Psen1 (Gene ID 19164) has 1 RefSeq proteins; selected NP_001349200.1 (467 aa, score 2243; best candidate score 2243).
- PSEN1: Rattus norvegicus ortholog Psen1 (Gene ID 29192) has 9 RefSeq proteins; selected NP_062036.2 (468 aa, score 2214; best candidate score 2214).
- PSEN1: Canis lupus familiaris ortholog PSEN1 (Gene ID 403408) has 6 RefSeq proteins; selected NP_001002943.2 (466 aa, score 2244; best candidate score 2244).
- PSEN1: Bos taurus ortholog PSEN1 (Gene ID 282705) has 10 RefSeq proteins; selected NP_777146.1 (478 aa, score 2185; best candidate score 2209).
- PSEN2: Mus musculus ortholog Psen2 (Gene ID 19165) has 6 RefSeq proteins; selected NP_001122077.1 (448 aa, score 2244; best candidate score 2244).
- PSEN2: Rattus norvegicus ortholog Psen2 (Gene ID 81751) has 45 RefSeq proteins; selected NP_112349.2 (448 aa, score 2227; best candidate score 2227).
- PSEN2: Canis lupus familiaris ortholog PSEN2 (Gene ID 490382) has 5 RefSeq proteins; selected NP_001271394.1 (448 aa, score 2230; best candidate score 2238).
- PSEN2: Bos taurus ortholog PSEN2 (Gene ID 282010) has 8 RefSeq proteins; selected NP_776865.2 (449 aa, score 2268; best candidate score 2268).
- FASTA validation passed: amino-acid alphabet only, no duplicate accessions, lengths within 0.5-2x of human, one sequence type (protein) per analysis set, headers carry accession and [species]. No cleaned copies were needed, so 02_raw_sequences/cleaned/ is not used.
- APP: BLASTp (BLASTP 2.17.0+, db refseq_protein) top hit per species checked against the ortholog gene (paralog check). Raw output APP_blastp_NP_000475.1_refseq_protein.json accessed 2026-10-29.
- PSEN1: BLASTp submitted, RID EKUTD9NR016.
