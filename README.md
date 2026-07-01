# VCF Trio Analysis — Structural Variant Screening & Annotation

Analysis of **structural variants (SVs)** called by [Delly](https://github.com/dellytools/delly)
in a parent–offspring **trio**, including a **de novo mutation screen**, sample-role
inference (father / mother / child), and **functional annotation** of candidate variants
with [ANNOVAR](https://annovar.openbioinformatics.org/).

Developed as a Next-Generation Sequencing (NGS) analysis project.

> **Data availability:** The raw VCF used in this project is access-restricted and is
> **not** included in this repository. See [`data/README.md`](data/README.md) for the exact
> input/output files the pipeline expects and how to reproduce the analysis with your own data.

---

## Contents

- [Background](#background)
- [Objective](#objective)
- [Repository structure](#repository-structure)
- [Pipeline overview](#pipeline-overview)
- [Key results](#key-results)
- [Reproducing the analysis](#reproducing-the-analysis)
- [Dependencies](#dependencies)
- [References](#references)
- [License](#license)

---

## Background

A **VCF (Variant Call Format)** file lists differences between sequenced samples and a
reference genome. Here the variants are **structural variants** (events larger than ~50 bp,
mainly insertions and deletions) produced by the Delly caller against the **hg38** human
reference.

Two ideas drive the analysis:

- **De novo mutations** are variants present in a child but in neither parent. They arise in
  the parental germline or early embryo and are a known cause of rare genetic disorders, so
  screening a trio for them is a standard clinical-genomics task.
- **Variant annotation** turns raw coordinates into biological meaning — which gene is hit,
  what the consequence is, how common the variant is in the population, and whether it has any
  reported clinical significance.

## Objective

> Analyse the variants in the Delly trio VCF, identify potential de novo mutations, and
> functionally annotate variants of interest to evaluate their biological importance.

---

## Repository structure

```
vcf-trio-analysis/
├── README.md                      # You are here
├── LICENSE                        # MIT license
├── CITATION.cff                   # How to cite this work (optional, GitHub-rendered)
├── .gitignore                     # Keeps raw data and generated artifacts out of version control
├── analysis/
│   └── vcf_trio_analysis.Rmd      # Full analysis: Bash preprocessing + R visualisation
├── data/
│   └── README.md                  # Expected data layout (no raw data committed)
├── results/
│   └── .gitkeep                   # Placeholder; generated figures/tables land here
└── docs/
    └── README.md                  # Where to put supporting documentation
```

---

## Pipeline overview

The workflow moves from raw calls to annotated, interpreted variants in five stages. The full,
runnable code lives in [`analysis/vcf_trio_analysis.Rmd`](analysis/vcf_trio_analysis.Rmd).

**1. Preprocess & subset the VCF (Bash / awk)**
Reduce the multi-sample Delly VCF to just the three trio samples (`HG00512`, `HG00513`,
`HG00514`), inspect the header, count variants, and split them by `FILTER` status (`PASS`
vs `LowQual`).

**2. Clean & screen for de novo mutations (awk)**
Remove records with ambiguous reference bases (`N`), then select sites where **both parents are
homozygous reference (`0/0`)** and the **child carries a non-reference, non-missing genotype**.

**3. Infer sample roles (awk + R)**
Build a genotype table to pick out the child, then use **X-chromosome read depth** to
distinguish the parents: the sample with roughly half the chrX read count carries a single X
and is inferred to be the father.

**4. Annotate variants (ANNOVAR, hg38)**
Convert the cleaned VCF to ANNOVAR `.avinput`, then annotate against four databases:

| Database          | Purpose                                                        |
|-------------------|----------------------------------------------------------------|
| `refGeneWithVer`  | Gene-based annotation (which gene, and the functional effect)   |
| `cytoBand`        | Chromosome band (Giemsa-stained cytogenetic location)          |
| `gnomad41_genome` | Population allele frequencies (gnomAD v4.1 whole genome)        |
| `clinvar_20250721`| Reported clinical significance (ClinVar)                       |

**5. Summarise & visualise (R)**
Load the ANNOVAR CSV in R and produce the summary tables and plots (variant type,
genomic-region distribution, per-chromosome counts, insertion/deletion split).

---

## Key results

- **Variant counts:** 17,312 total variants — **17,152 PASS**, 159 LowQual. The set is
  dominated by structural indels: **9,444 deletions** and **5,315 insertions**; no simple
  bi-allelic SNPs and no multi-allelic sites were detected.
- **De novo screen:** the three candidate sites initially found all contained ambiguous (`N`)
  bases. After filtering those out, **no de novo mutations remained**.
- **Sample roles:** genotype comparison identified **HG00514 as the child**; the lower
  X-chromosome read depth identified **HG00512 as the father** (HG00513 and the child inferred
  female).
- **Genomic distribution:** most variants are **intergenic (~54.6%)** or **intronic (~35.8%)**,
  with only a small exonic fraction.
- **Variants of interest:** eight exonic frameshift / non-frameshift indels across chromosomes
  1, 2, 3, 4, 9, 14 and 18, in the genes *SEC16B*, *SULT1C3*, *BARD1*, *MUC4*, *RNF212*,
  *PRUNE2*, *OR4L1* and *EMILIN2*. Gene-level intolerance metrics (pLI, RVIS, GDI) were used to
  gauge how sensitive each gene is to variation.
- **Highlight — *EMILIN2*:** a rare frameshift variant (`p.K787Sfs*126`, population allele
  frequency < 0.1%) predicted to cause loss of function in an extracellular-matrix glycoprotein.
- **Conclusion:** no clinically relevant variants were identified. *EMILIN2* carries a rare
  variant, but there is no evidence that any of the candidate genes have a harmful effect in
  this trio.

---

## Reproducing the analysis

You will need your own trio VCF (see [`data/README.md`](data/README.md) for the expected
layout). At a high level:

1. Place the raw Delly VCF where the notebook expects it (see `data/README.md`).
2. Run the Bash/awk preprocessing and de novo steps (embedded in the `.Rmd`, or lift them into
   a shell script).
3. Download the ANNOVAR databases listed above and run the annotation to produce the CSV.
4. Knit the notebook to generate the tables and figures:

   ```r
   rmarkdown::render("analysis/vcf_trio_analysis.Rmd")
   ```

> Some Bash chunks in the notebook are commented out on purpose — they document
> environment-specific steps (fetching the source data, running ANNOVAR) that won't execute in
> an arbitrary setup. Treat them as a reproducible record of what was run.

---

## Dependencies

**Command-line**
- A POSIX shell with `awk`, `grep`, `sed`, `sort`, `uniq`, `wc`
- [ANNOVAR](https://annovar.openbioinformatics.org/) (with the hg38 databases above)

**R (≥ 4.x)** — install with:

```r
install.packages(c(
  "tidyverse", "dplyr", "ggplot2",
  "RColorBrewer", "viridis", "stringr", "rmarkdown"
))
```

---

## References

- Johnston, M.O. (2025). *Mutations and New Variation: Overview.* eLS.
  https://doi.org/10.1002/9780470015902.a0029422
- Nicolas, G. & Veltman, J.A. (2018). *The role of de novo mutations in adult-onset
  neurodegenerative disorders.* Acta Neuropathologica, 137(2), 183–207.
  https://doi.org/10.1007/s00401-018-1939-3
- Wang, W., Corominas, R. & Lin, G.N. (2019). *De novo Mutations From Whole Exome Sequencing
  in Neurodevelopmental and Psychiatric Disorders.* Frontiers in Genetics, 10.
  https://doi.org/10.3389/fgene.2019.00258
- Petrovski, S. et al. (2013). *Genic intolerance to functional variation and the
  interpretation of personal genomes* (RVIS). PLoS Genetics, 9(8), e1003709.
  https://doi.org/10.1371/journal.pgen.1003709
- Alyousfi, D. et al. (2018). *Gene-specific metrics to inform gene prioritisation.* Briefings
  in Functional Genomics, 18(1), 23–29. https://doi.org/10.1093/bfgp/ely033
- ANNOVAR — https://annovar.openbioinformatics.org/
- Delly — https://github.com/dellytools/delly

---

## License

Released under the [MIT License](LICENSE). See the license file for details.
