### Rodent Research-1 CASIS experiment, microgravity-associated muscle wasting.
- OSD-47 Mouse liver transcriptomic, proteomic, epigenomic and histology data

[Data source](https://osdr.nasa.gov/bio/repo/data/studies/OSD-47)

- **FLT**: Dissected on orbit 21/22 days after launch
- **GC**: Age-matched Ground Controls
- **BSL**: Basal controls (euthanized at time of launch)

## Snakemake Pipeline Bulk RNA-seq

```mermaid
flowchart TD
	classDef control fill:#F3F4F6,stroke:#4B5563,color:#111827,stroke-width:1.2px;
	classDef config fill:#FFF4E6,stroke:#D17B0F,color:#111827,stroke-width:1.2px;
	classDef input fill:#EEF4FF,stroke:#4C78A8,color:#111827,stroke-width:1.1px;
	classDef stage fill:#F9F6E7,stroke:#B8871B,color:#111827,stroke-width:1.2px;
	classDef output fill:#EDF7ED,stroke:#4E8F5C,color:#111827,stroke-width:1.1px;

	snakefile[Snakefile<br/>rule all]:::control
	config[config/config.yaml<br/>pipeline settings]:::config
	input[GeneLab RNA-seq archives<br/>HISAT2 index + Ensembl GTF]:::input

	archive[Archive + manifests<br/>extract and build sample metadata]:::stage
	trim[Trim + read QC<br/>fastp, FastQC, MultiQC]:::stage
	align[Merge + align<br/>HISAT2, samtools, BAM QC]:::stage
	counts[Count genes<br/>featureCounts matrix]:::stage
	de[Normalize + contrasts<br/>PyDESeq2 outputs]:::stage
	report[Plots + report<br/>qc_diff_exp outputs]:::output

	config -.-> snakefile
	input --> archive --> trim --> align --> counts --> de --> report
	snakefile --> archive
```

### Example volcano plots from normalized counts (Bulk RNA-seq):

![FLT_BASE](plots/condition_FLT_vs_BSL.volcano.png)
![FLT_GC](plots/condition_FLT_vs_GC.volcano.png)

## Proteomics

1. `.raw` $\rightarrow$ `mzML`
2. `mzML` + [FragPipe](https://fragpipe.nesvilab.org/) $\rightarrow$ quantitative proteomics analyses at different resolutions (e.g., gene-leve, protein-level)

Sample-level tissue-marker QC scores (mean per-gene z-scores across muscle and liver marker panels):

| sample | liver | muscle | muscle_minus_liver |
| --- | ---: | ---: | ---: |
| BSL_Rep1 | -0.843 | -0.636 | 0.208 |
| BSL_Rep2 | -0.739 | -0.633 | 0.106 |
| BSL_Rep3 | 0.079 | -0.437 | -0.516 |
| FLT_Rep1 | 0.852 | -0.572 | -1.423 |
| FLT_Rep2 | -0.128 | 0.034 | 0.162 |
| FLT_Rep3 | 0.369 | -0.431 | -0.800 |
| GC_Rep1 | -0.301 | 2.371 | 2.672 |
| GC_Rep2 | 0.439 | -0.350 | -0.789 |
| GC_Rep3 | 0.273 | 0.654 | 0.380 |

* GC_Rep1 possibly contaminated with muscle tissue
* Exclude all GC_rep for this run

![heatmap](plots/heatmap_proteome_no_GC.png)


