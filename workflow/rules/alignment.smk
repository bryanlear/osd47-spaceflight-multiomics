TRIMMING = config.get("trimming", {})
TRIMMING_ENABLED = TRIMMING.get("enabled", False)
ALIGNMENT = config.get("alignment", {})
ALIGNMENT_ENABLED = ALIGNMENT.get("enabled", False)
ALIGNMENT_TOOL = ALIGNMENT.get("tool", "hisat2")
ARCHIVE_IDS = [archive["id"] for archive in config["archives"]]


if ALIGNMENT_ENABLED and not TRIMMING_ENABLED:
	raise ValueError("Alignment requires trimming.enabled: true because it consumes results/trimmed outputs.")

if ALIGNMENT_ENABLED and ALIGNMENT_TOOL != "hisat2":
	raise ValueError(f"Unsupported alignment tool: {ALIGNMENT_TOOL}")

if ALIGNMENT_ENABLED and not ALIGNMENT.get("hisat2_index_base"):
	raise ValueError("Alignment requires alignment.hisat2_index_base to be set.")


HISAT2_INDEX_BASE = ALIGNMENT.get("hisat2_index_base")
HISAT2_INDEX_FILES = (
	[f"{HISAT2_INDEX_BASE}.{index}.ht2" for index in range(1, 9)]
	if ALIGNMENT_ENABLED and HISAT2_INDEX_BASE
	else []
)


rule merge_trimmed_fastq:
	input:
		manifest=config["outputs"]["lane_manifest"],
		trimmed_marker="results/trimmed/.trimmed.ok"
	output:
		merged="results/merged_trimmed/{archive_id}.trimmed.fastq.gz"
	params:
		trimmed_root="results/trimmed"
	shell:
		r"""
		set -euo pipefail
		mkdir -p "$(dirname "{output.merged}")"
		inputs=()
		while IFS=$'\t' read -r archive_id condition replicate lane_fastq lane_name; do
			if [[ "${{archive_id}}" != "{wildcards.archive_id}" ]]; then
				continue
			fi

			lane_base="${{lane_name%.fastq.gz}}"
			if [[ "${{lane_base}}" == "${{lane_name}}" ]]; then
				lane_base="${{lane_name%.fq.gz}}"
			fi

			inputs+=("{params.trimmed_root}/${{archive_id}}/${{lane_base}}.trimmed.fastq.gz")
		done < <(tail -n +2 "{input.manifest}")

		if [[ "${{#inputs[@]}}" -eq 0 ]]; then
			echo "No trimmed lanes found for {wildcards.archive_id}" >&2
			exit 1
		fi

		cat "${{inputs[@]}}" > "{output.merged}"
		"""


rule align_trimmed_sample:
	input:
		merged="results/merged_trimmed/{archive_id}.trimmed.fastq.gz",
		index=HISAT2_INDEX_FILES
	output:
		bam="results/hisat2_alignments/{archive_id}.sorted.bam",
		bai="results/hisat2_alignments/{archive_id}.sorted.bam.bai",
		summary="results/hisat2_alignments/{archive_id}.hisat2.summary.txt",
		stderr="results/hisat2_alignments/{archive_id}.hisat2.stderr.log",
		flagstat="results/hisat2_alignments/{archive_id}.flagstat.txt"
	threads:
		lambda wildcards: ALIGNMENT.get("threads", 8)
	params:
		index_base=HISAT2_INDEX_BASE,
		rna_strandness=lambda wildcards: ALIGNMENT.get("rna_strandness", "R"),
		sort_threads=lambda wildcards: ALIGNMENT.get("sort_threads", 4)
	shell:
		r"""
		set -euo pipefail
		mkdir -p "$(dirname "{output.bam}")"
		hisat2 -p {threads} \
		  --rna-strandness {params.rna_strandness} \
		  -x "{params.index_base}" \
		  -U "{input.merged}" \
		  --summary-file "{output.summary}" \
		  2> "{output.stderr}" \
		| samtools sort -@ {params.sort_threads} -o "{output.bam}"

		samtools index "{output.bam}" "{output.bai}"
		samtools flagstat "{output.bam}" > "{output.flagstat}"
		"""