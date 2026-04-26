configfile: "config/config.yaml"

TRIMMING = config.get("trimming", {})
TRIMMING_ENABLED = TRIMMING.get("enabled", False)
STUDY_SLUG = str(config["study"]["accession"]).lower().replace("-", "")
TRIMMED_MULTIQC_HTML = f"results/qc/trimmed_multiqc/{STUDY_SLUG}_trimmed_multiqc.html"

BASE_TARGETS = [
    config["outputs"]["archive_manifest"],
    config["outputs"]["lane_manifest"],
    *expand(
        "results/extracted/{archive}/.extracted.ok",
        archive=[archive["id"] for archive in config["archives"]],
    ),
]
TRIMMING_TARGETS = (
    ["results/trimmed/.trimmed.ok", TRIMMED_MULTIQC_HTML] if TRIMMING_ENABLED else []
)

include: "workflow/rules/archive.smk"
include: "workflow/rules/trimming.smk"


rule all:
    input:
        BASE_TARGETS + TRIMMING_TARGETS