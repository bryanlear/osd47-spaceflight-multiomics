### Rodent Research-1 CASIS experiment, microgravity-associated muscle wasting.
- OSD-47 Mouse liver transcriptomic, proteomic, epigenomic and histology data

[Data source](https://osdr.nasa.gov/bio/repo/data/studies/OSD-47)

- **FLT**: Dissected on orbit 21/22 days after launch
- **GC**: Age-matched Ground Controls
- **BC**: Basal controls (euthanized at time of launch)

### Snakemake Pipeline

```mermaid
flowchart TD
	classDef control fill:#F3F4F6,stroke:#4B5563,color:#111827,stroke-width:1.2px;
	classDef config fill:#FFF4E6,stroke:#D17B0F,color:#111827,stroke-width:1.2px;
	classDef input fill:#EEF4FF,stroke:#4C78A8,color:#111827,stroke-width:1.1px;
	classDef stage fill:#F9F6E7,stroke:#B8871B,color:#111827,stroke-width:1.2px;
	classDef output fill:#EDF7ED,stroke:#4E8F5C,color:#111827,stroke-width:1.1px;

	snakefile[Snakefile<br/>rule all collects final targets]:::control
	config[config/config.yaml<br/>toggles, tool settings, output paths]:::config

	raw_archives[GeneLab RNA-seq archives<br/>GLDS-47 tar.gz inputs]:::input
	references[Reference inputs<br/>HISAT2 index and Ensembl GTF]:::input

	archive_stage[Archive stage<br/>archive.smk<br/>build_archive_manifest<br/>extract_archive<br/>build_lane_manifest]:::stage
	archive_outputs[Metadata and extraction outputs<br/>archive_manifest.tsv<br/>lane_manifest.tsv<br/>results/extracted/*/.extracted.ok]:::output

	trim_stage[Trimming stage<br/>trimming.smk<br/>trim_lanes<br/>run_trimmed_fastqc<br/>run_trimmed_multiqc]:::stage
	trimmed_reads[Trimmed reads<br/>results/trimmed/archive_id/lane.trimmed.fastq.gz]:::output
	trimmed_qc[Trimmed-read QC<br/>FastQC outputs and MultiQC HTML]:::output

	align_stage[Alignment stage<br/>alignment.smk<br/>merge_trimmed_fastq<br/>align_trimmed_sample]:::stage
	alignment_outputs[Merged FASTQ and alignment outputs<br/>results/merged_trimmed/*<br/>results/hisat2_alignments/*.sorted.bam<br/>plus BAI, HISAT2 summary, flagstat]:::output

	count_stage[Counting stage<br/>alignment.smk<br/>run_featurecounts]:::stage
	count_outputs[Gene count matrix<br/>results/counts/osd47_featureCounts.txt<br/>and .summary]:::output

	de_stage[Differential-expression stage<br/>differential_expression.smk<br/>run_pydeseq2]:::stage
	de_outputs[Normalized DE outputs<br/>results/pydeseq2_all/<br/>sample_metadata.tsv<br/>counts_matrix.samples_x_genes.tsv<br/>contrast tables]:::output

	qc_stage[QC and reporting stage<br/>differential_expression.smk<br/>run_qc_diff_exp_analysis]:::stage
	qc_outputs[Final plots and report<br/>results/pydeseq2_all/qc_diff_exp/<br/>library QC, PCA, heatmap, volcano plots<br/>top_up and top_down tables, report.md]:::output

	config -.-> snakefile
	config -.-> trim_stage
	config -.-> align_stage
	config -.-> count_stage
	config -.-> de_stage
	config -.-> qc_stage

	raw_archives --> archive_stage --> archive_outputs
	archive_outputs --> trim_stage --> trimmed_reads
	trim_stage --> trimmed_qc
	archive_outputs --> align_stage
	trimmed_reads --> align_stage
	references --> align_stage
	align_stage --> alignment_outputs --> count_stage --> count_outputs
	references --> count_stage
	count_outputs --> de_stage --> de_outputs --> qc_stage --> qc_outputs
	references --> qc_stage

	snakefile --> archive_stage
	snakefile --> trim_stage
	snakefile --> align_stage
	snakefile --> count_stage
	snakefile --> de_stage
	snakefile --> qc_stage
```

### Examples volcano plots from normalized data

![FLT_BASE](plots/condition_FLT_vs_BSL.volcano.png)
![FLT_GC](plots/condition_FLT_vs_GC.volcano.png)