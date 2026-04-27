configfile: "config/config.yaml"

TRIMMING = config.get("trimming", {})
TRIMMING_ENABLED = TRIMMING.get("enabled", False)
ALIGNMENT = config.get("alignment", {})
ALIGNMENT_ENABLED = ALIGNMENT.get("enabled", False)
ARCHIVE_IDS = [archive["id"] for archive in config["archives"]]
STUDY_SLUG = str(config["study"]["accession"]).lower().replace("-", "")
TRIMMED_MULTIQC_HTML = f"results/qc/trimmed_multiqc/{STUDY_SLUG}_trimmed_multiqc.html"

BASE_TARGETS = [
    config["outputs"]["archive_manifest"],
    config["outputs"]["lane_manifest"],
    *expand(
        "results/extracted/{archive}/.extracted.ok",
        archive=ARCHIVE_IDS,
    ),
]
TRIMMING_TARGETS = (
    ["results/trimmed/.trimmed.ok", TRIMMED_MULTIQC_HTML] if TRIMMING_ENABLED else []
)
ALIGNMENT_TARGETS = (
    expand("results/hisat2_alignments/{archive}.sorted.bam", archive=ARCHIVE_IDS)
    if ALIGNMENT_ENABLED
    else []
)

include: "workflow/rules/archive.smk"
include: "workflow/rules/trimming.smk"
include: "workflow/rules/alignment.smk"


rule all:
    input:
        BASE_TARGETS + TRIMMING_TARGETS + ALIGNMENT_TARGETS