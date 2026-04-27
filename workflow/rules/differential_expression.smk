FEATURE_COUNTS = config.get("feature_counts", {})
FEATURE_COUNTS_ENABLED = FEATURE_COUNTS.get("enabled", False)
PYDESEQ2 = config.get("pydeseq2", {})
PYDESEQ2_ENABLED = PYDESEQ2.get("enabled", False)
CONDITIONS = list(dict.fromkeys(str(archive["condition"]) for archive in config["archives"]))


def build_pydeseq2_contrast_pairs():
	design_factor = PYDESEQ2.get("design_factor", "condition")
	if design_factor != "condition":
		raise ValueError(
			"Snakemake PyDESeq2 integration currently supports pydeseq2.design_factor: condition only."
		)

	reference_level = str(PYDESEQ2.get("reference_level", "BSL"))
	contrast_level = str(PYDESEQ2.get("contrast_level", "FLT"))

	if reference_level not in CONDITIONS:
		raise ValueError(
			f"PyDESeq2 reference level '{reference_level}' is not present in archive conditions: {CONDITIONS}"
		)
	if contrast_level not in CONDITIONS:
		raise ValueError(
			f"PyDESeq2 contrast level '{contrast_level}' is not present in archive conditions: {CONDITIONS}"
		)

	ordered_levels = [reference_level] + [
		level for level in CONDITIONS if level != reference_level
	]

	if not PYDESEQ2.get("all_pairwise", True):
		return [(contrast_level, reference_level)]

	contrast_pairs = [
		(contrast_level, level) for level in ordered_levels if level != contrast_level
	]
	remaining_levels = [level for level in ordered_levels if level != contrast_level]
	for lower_index, left_level in enumerate(remaining_levels[:-1]):
		for right_level in remaining_levels[lower_index + 1 :]:
			contrast_pairs.append((right_level, left_level))
	return contrast_pairs


if PYDESEQ2_ENABLED and not FEATURE_COUNTS_ENABLED:
	raise ValueError("pydeseq2 requires feature_counts.enabled: true because it consumes the count matrix.")


FEATURE_COUNTS_OUTPUT = FEATURE_COUNTS.get(
	"output",
	f"results/counts/{config['study']['accession'].lower().replace('-', '')}_featureCounts.txt",
)
PYDESEQ2_OUTPUT_DIR = PYDESEQ2.get("output_dir", "results/pydeseq2")
PYDESEQ2_CONTRAST_STEMS = (
	[
		f"{PYDESEQ2.get('design_factor', 'condition')}_{contrast}_vs_{reference}"
		for contrast, reference in build_pydeseq2_contrast_pairs()
	]
	if PYDESEQ2_ENABLED
	else []
)
PYDESEQ2_TARGETS = (
	[
		f"{PYDESEQ2_OUTPUT_DIR}/sample_metadata.tsv",
		f"{PYDESEQ2_OUTPUT_DIR}/counts_matrix.samples_x_genes.tsv",
		f"{PYDESEQ2_OUTPUT_DIR}/contrast_manifest.tsv",
		*[f"{PYDESEQ2_OUTPUT_DIR}/{stem}.tsv" for stem in PYDESEQ2_CONTRAST_STEMS],
		*[
			f"{PYDESEQ2_OUTPUT_DIR}/{stem}.significant.tsv"
			for stem in PYDESEQ2_CONTRAST_STEMS
		],
	]
	if PYDESEQ2_ENABLED
	else []
)


rule run_pydeseq2:
	input:
		counts=FEATURE_COUNTS_OUTPUT,
		counts_summary=f"{FEATURE_COUNTS_OUTPUT}.summary"
	output:
		metadata=f"{PYDESEQ2_OUTPUT_DIR}/sample_metadata.tsv",
		counts_matrix=f"{PYDESEQ2_OUTPUT_DIR}/counts_matrix.samples_x_genes.tsv",
		contrast_manifest=f"{PYDESEQ2_OUTPUT_DIR}/contrast_manifest.tsv",
		contrasts=expand(f"{PYDESEQ2_OUTPUT_DIR}/{{stem}}.tsv", stem=PYDESEQ2_CONTRAST_STEMS),
		significant=expand(
			f"{PYDESEQ2_OUTPUT_DIR}/{{stem}}.significant.tsv",
			stem=PYDESEQ2_CONTRAST_STEMS,
		)
	threads:
		lambda wildcards: PYDESEQ2.get("n_cpus", 8)
	params:
		config_path="config/config.yaml",
		output_dir=PYDESEQ2_OUTPUT_DIR,
		design_factor=lambda wildcards: PYDESEQ2.get("design_factor", "condition"),
		reference_level=lambda wildcards: PYDESEQ2.get("reference_level", "BSL"),
		contrast_level=lambda wildcards: PYDESEQ2.get("contrast_level", "FLT"),
		all_pairwise_flag="--all-pairwise" if PYDESEQ2.get("all_pairwise", True) else "",
		min_count=lambda wildcards: PYDESEQ2.get("min_count", 10),
		min_samples=lambda wildcards: PYDESEQ2.get("min_samples", 3),
		alpha=lambda wildcards: PYDESEQ2.get("alpha", 0.05)
	shell:
		r"""
		set -euo pipefail
		mkdir -p "{params.output_dir}"
		python workflow/scripts/run_pydeseq2.py \
		  --counts "{input.counts}" \
		  --config "{params.config_path}" \
		  --output-dir "{params.output_dir}" \
		  --design-factor "{params.design_factor}" \
		  --reference-level "{params.reference_level}" \
		  --contrast-level "{params.contrast_level}" \
		  --min-count {params.min_count} \
		  --min-samples {params.min_samples} \
		  --alpha {params.alpha} \
		  --n-cpus {threads} \
		  {params.all_pairwise_flag}
		"""