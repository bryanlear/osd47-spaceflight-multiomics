import shlex


TRIMMING = config.get("trimming", {})
TRIMMING_ENABLED = TRIMMING.get("enabled", False)
TRIMMING_TOOL = TRIMMING.get("tool", "fastp")
STUDY_SLUG = str(config["study"]["accession"]).lower().replace("-", "")
TRIMMED_MULTIQC_STEM = f"{STUDY_SLUG}_trimmed_multiqc"


if TRIMMING_ENABLED and TRIMMING_TOOL != "fastp":
	raise ValueError(f"Unsupported trimming tool: {TRIMMING_TOOL}")


def fastp_extra_args():
	args = []

	adapter_sequence = TRIMMING.get("adapter_sequence")
	if adapter_sequence:
		args.extend(["-a", adapter_sequence])

	if TRIMMING.get("cut_right", False):
		args.extend(
			[
				"-r",
				"--cut_right_window_size",
				str(TRIMMING.get("cut_window_size", 4)),
				"--cut_right_mean_quality",
				str(TRIMMING.get("cut_mean_quality", 20)),
			]
		)

	if TRIMMING.get("poly_x", False):
		args.extend(["-x", "--poly_x_min_len", str(TRIMMING.get("poly_x_min_len", 10))])

	return " ".join(shlex.quote(arg) for arg in args)


rule trim_lanes:
	input:
		manifest=config["outputs"]["lane_manifest"]
	output:
		marker="results/trimmed/.trimmed.ok"
	threads:
		lambda wildcards: TRIMMING.get("threads", 4)
	params:
		trimmed_root="results/trimmed",
		report_root="results/qc/trimmed_fastp",
		extra_args=fastp_extra_args(),
		qualified_quality=lambda wildcards: TRIMMING.get("qualified_quality_phred", 20),
		unqualified_percent=lambda wildcards: TRIMMING.get("unqualified_percent_limit", 30),
		min_length=lambda wildcards: TRIMMING.get("min_length", 20)
	shell:
		r"""
		set -euo pipefail
		mkdir -p {params.trimmed_root} {params.report_root}
		tail -n +2 "{input.manifest}" | while IFS=$'\t' read -r archive_id condition replicate lane_fastq lane_name; do
			lane_base="${{lane_name%.fastq.gz}}"
			if [[ "${{lane_base}}" == "${{lane_name}}" ]]; then
				lane_base="${{lane_name%.fq.gz}}"
			fi

			mkdir -p "{params.trimmed_root}/${{archive_id}}" "{params.report_root}/${{archive_id}}"

			fastp \
			  -i "${{lane_fastq}}" \
			  -o "{params.trimmed_root}/${{archive_id}}/${{lane_base}}.trimmed.fastq.gz" \
			  -q {params.qualified_quality} \
			  -u {params.unqualified_percent} \
			  -l {params.min_length} \
			  -h "{params.report_root}/${{archive_id}}/${{lane_base}}.fastp.html" \
			  -j "{params.report_root}/${{archive_id}}/${{lane_base}}.fastp.json" \
			  -w {threads} {params.extra_args}
		done
		touch "{output.marker}"
		"""


rule run_trimmed_fastqc:
	input:
		manifest=config["outputs"]["lane_manifest"],
		trimmed_marker="results/trimmed/.trimmed.ok"
	output:
		marker="results/qc/trimmed_fastqc/.fastqc.ok"
	threads:
		lambda wildcards: TRIMMING.get("threads", 4)
	params:
		trimmed_root="results/trimmed",
		report_root="results/qc/trimmed_fastqc"
	shell:
		r"""
		set -euo pipefail
		mkdir -p {params.report_root}
		rm -f "{output.marker}"
		tail -n +2 "{input.manifest}" | while IFS=$'\t' read -r archive_id condition replicate lane_fastq lane_name; do
			lane_base="${{lane_name%.fastq.gz}}"
			if [[ "${{lane_base}}" == "${{lane_name}}" ]]; then
				lane_base="${{lane_name%.fq.gz}}"
			fi

			fastqc \
			  --threads {threads} \
			  --outdir "{params.report_root}" \
			  "{params.trimmed_root}/${{archive_id}}/${{lane_base}}.trimmed.fastq.gz"
		done
		touch "{output.marker}"
		"""


rule run_trimmed_multiqc:
	input:
		trimmed_marker="results/trimmed/.trimmed.ok",
		fastqc_marker="results/qc/trimmed_fastqc/.fastqc.ok"
	output:
		html=f"results/qc/trimmed_multiqc/{TRIMMED_MULTIQC_STEM}.html",
		data_dir=directory(f"results/qc/trimmed_multiqc/{TRIMMED_MULTIQC_STEM}_data")
	params:
		report_root="results/qc/trimmed_multiqc",
		report_name=f"{TRIMMED_MULTIQC_STEM}.html"
	shell:
		r"""
		set -euo pipefail
		mkdir -p {params.report_root}
		multiqc \
		  --force \
		  --outdir "{params.report_root}" \
		  --filename "{params.report_name}" \
		  "results/qc/trimmed_fastqc" \
		  "results/qc/trimmed_fastp"
		"""
