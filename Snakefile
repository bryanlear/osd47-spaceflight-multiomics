configfile: "config/config.yaml"

TRIMMING = config.get("trimming", {})
TRIMMING_ENABLED = TRIMMING.get("enabled", False)

BASE_TARGETS = [
    config["outputs"]["archive_manifest"],
    config["outputs"]["lane_manifest"],
    *expand(
        "results/extracted/{archive}/.extracted.ok",
        archive=[archive["id"] for archive in config["archives"]],
    ),
]
TRIMMING_TARGETS = ["results/trimmed/.trimmed.ok"] if TRIMMING_ENABLED else []

include: "workflow/rules/archive.smk"
include: "workflow/rules/trimming.smk"


rule all:
    input:
        BASE_TARGETS + TRIMMING_TARGETS