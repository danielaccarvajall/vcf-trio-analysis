# Data

## Data availability statement

The analysis is based on a multi-sample **Delly** structural-variant VCF (hg38) covering a
parent–offspring trio. The samples analysed are `HG00512`, `HG00513`, and `HG00514`. The raw
VCF is access-restricted and cannot be redistributed here; anyone with access to an equivalent
trio VCF can reproduce the workflow using the code in
[`../analysis/vcf_trio_analysis.Rmd`](../analysis/vcf_trio_analysis.Rmd).

## Expected files

Place your input here (kept untracked by git) and the pipeline will generate the rest:

| File                       | Stage       | Description                                                        |
|----------------------------|-------------|--------------------------------------------------------------------|
| `DellyVariation.vcf`       | input       | Raw multi-sample Delly SV calls (hg38). **Provide this yourself.** |
| `DellyTrio1.vcf`           | generated   | VCF subset to the three trio samples                               |
| `CleanDellyTrio1Pass.vcf`  | generated   | Trio VCF with ambiguous-base (`N`) records removed                 |
| `genotypes.txt`            | generated   | Per-sample genotype table (PASS variants)                          |
| `ReadCounts.txt` / `CleanReadCounts.txt` | generated | Per-sample read counts, used for parent inference    |
| `CleanedDT1Pass.avinput`   | generated   | ANNOVAR input converted from the cleaned VCF                       |
| `annovar_output.csv`       | generated   | ANNOVAR annotation output loaded into R                            |

## Minimal VCF reference (for anyone new to the format)

A VCF has a header (lines starting with `#`) followed by one variant per line. The core columns:

```
CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  <sample columns...>
```

- **REF / ALT** — reference and alternate alleles (here, the alleles flank structural
  insertions/deletions).
- **FILTER** — quality flag; this project keeps `PASS` and separates out `LowQual`.
- **FORMAT + sample columns** — per-sample fields such as the genotype `GT` (`0/0`, `0/1`,
  `1/1`, or `./.` for missing) and read-count fields used here to infer sample sex/role.
- **INFO** — variant-level annotations; Delly records the SV type here (e.g. `SVTYPE=DEL`,
  `SVTYPE=INS`).
