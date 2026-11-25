# Validation of the additions in Phagree 1.1

The additions of version 1.1 do not change the outputs of the study (see "Regression" below). They were
checked on the study cohort with the scripts in `scripts/validation/`; the numbers below are in
`results/validation/`. Local PHROGs hits were available for 946,140 of the 1,333,849 distinct proteins of the
cohort (the rest were searched on another machine); function checks use these proteins.

## Search space of PHROGs E-values

HMMER computes E-values for a search space that, by default, is the number of proteins searched together
(Z) and, for domains, the number of significant proteins per profile (domZ). In the study, proteins were
searched in blocks of 20,000; in Phagree 1.0.0 `annotate` searched all proteins of a run as one block, so the
functions of borderline proteins depended on how many genomes were annotated together (6 of 3,215 proteins
of 28 laboratory genomes).

From 1.1, Z = 20,000 and domZ = 13 (the hit-weighted median number of significant proteins per profile in
the study blocks).

- 300 proteins of a full study block searched alone: E-values and bitscores identical to the study for all
  783 hits.
- Functions agree with the study for 99.8% of 2,045 proteins (the study itself used a domZ that varied
  between blocks, so exact agreement is not possible).
- 60 proteins searched alone and inside a block of 2,000: identical hits; phage lambda annotated alone and
  together with another genome: identical results.
- Recomputing the study rule from the cached hits reproduces the published function table for 944,876 of
  944,883 proteins. The 7 exceptions lie exactly on a threshold (for example E = 1.0 × 10^-5 or profile
  coverage 0.496) and were classified in the published table from a search in another block, i.e. they are
  the block-size effect described above.

## Long-profile rule (`--function-rule extended`)

Some PHROGs profiles are much longer than their members, so a complete protein can cover less than half of
the profile; for example the large terminase subunit of phage lambda (gpA, 641 aa) covers 46% of
`phrog_15` (1,330 positions) and has no specific function under the study rule. With the extended rule, a
protein without a specific function under the study rule is classified again with hits that cover at least
80% of the protein and at least 30% of the profile; the study result is never changed when it is a specific
function or ambiguous (0 study functions changed in the cohort).

Seven threshold variants were compared by the share of named GenBank products that share a function word
with the PHROGs annotation (a coarse lexical check, not accuracy):

| Variant (protein / profile coverage) | Proteins rescued | Agreement | Agreement with synonyms |
|---|---|---|---|
| Study rule (baseline) | — | 0.796 | 0.841 |
| ≥ 0.9 / any | 52,428 | 0.721 | 0.767 |
| ≥ 0.95 / any | 49,151 | 0.718 | 0.765 |
| ≥ 0.9 / ≥ 0.2 | 45,722 | 0.728 | 0.774 |
| ≥ 0.9 / ≥ 0.3 | 37,223 | 0.737 | 0.788 |
| ≥ 0.95 / ≥ 0.3 | 34,974 | 0.733 | 0.785 |
| **≥ 0.8 / ≥ 0.3 (chosen)** | **39,262** | **0.740** | **0.790** |
| ≥ 0.95 / ≥ 0.4 | 19,683 | 0.749 | 0.772 |

The criterion fixed before the comparison (agreement ≥ 0.78 without synonyms) was not met by any variant.
In 60 random examples most disagreements were synonyms (peptidase M15 / endolysin, major capsid / major
head protein, dUTP diphosphatase / dUTPase, baseplate / tail protein) or uninformative GenBank names
(ORF013); with the same synonym list applied to both groups, rescued functions still agree 5 points less
often than study functions. The extended rule is therefore **not the default**: it is available as an
option, uses the variant with the highest agreement and also more rescued proteins than its closest
alternative (≥ 0.8 / ≥ 0.3), and every function found this way is marked in the outputs ("assigned with the
long-profile rule").

## Lysogeny markers

Labels: INPHARED2 `phageleads_lifestyle` (PhageLeads predictions, not experimental data); 197 genomes with
label "unknown" were excluded. Relatedness blocks were split in two halves by a hash of the block name:
development (8,478 genomes, 2,253 blocks) and validation (9,187 genomes, 2,180 blocks). A genome has a
marker when any protein (GenBank or caller model) has a specific PHROGs function of the marker group.

The first selection rule fixed in advance (groups present in ≤ 1% of lytic and ≥ 1% of temperate genomes
of the development half) excluded integrase (2.2% of lytic genomes) and selected CI-like repressor,
Arc-like repressor and recombinase; their validation sensitivity was 11.7%. It was replaced, using the
development half only, by greedy selection: add the marker group that most increases Youden's J
(sensitivity + specificity − 1) and stop when the gain is below 0.005. This selected integrase,
transcriptional repressor and anti-repressor. Validation half, 95% cluster-bootstrap CIs over blocks:

| Evidence | Sensitivity (temperate) | Specificity (lytic) |
|---|---|---|
| `strong`: integrase | 0.68 (0.58–0.78) | 0.981 (0.974–0.987) |
| `strong` or `weak`: integrase, transcriptional repressor or anti-repressor | 0.86 (0.81–0.91) | 0.91 (0.85–0.95) |

`none` means that no marker was found, not that the phage is virulent. Per-marker rates are in
`results/validation/lysogeny_marker_validation.json`.

## NCBI submission files

The feature tables and sequences of phage lambda (as a GenBank record and as a FASTA genome) pass
`table2asn` 1.29.324 (`-V vb`) without annotation errors or warnings; the only messages concern the
submission itself (no publication or submitter template, taxon ID not looked up offline), which the
submitter adds.

## Regression

With default options, the re-annotation of the 17,859 study records by Phagree 1.1 is byte-identical to
`results/reannotation/` of 1.0.0 for all five files (`records_reannotation.csv.gz`, `gene_proposals.csv.gz`,
`consensus_annotation.gff3.gz`, `sequence_quality.csv.gz`, `reannotation_summary.json`); the new files
(`consensus_genes.csv.gz`, `protein_function_extended.parquet`) are added next to them.

`phagree annotate` was run on seven genomes (phages lambda, phiX174, T4 and T7 and three laboratory phages,
MN510331.1, PQ474795.1 and OM960734.1 with 1,231 ambiguous bases), each as a GenBank record and as a FASTA
genome, with both function rules. All runs completed; lysogeny evidence was `strong` for lambda (integrase)
and `none` for T4, T7 and phiX174; the extended rule added 83 functions over the 14 genomes. The NCBI files
of all 14 genomes pass table2asn without annotation errors or warnings (phiX174, a circular genome, once its
viral lineage is given in the source information, as in a real submission).
