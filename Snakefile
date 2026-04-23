configfile: "config/config.yaml"

include: "workflow/rules/archive.smk"


rule all:
    input:
        config["outputs"]["archive_manifest"],
        config["outputs"]["lane_manifest"],
        expand(
            "results/extracted/{archive}/.extracted.ok",
            archive=[archive["id"] for archive in config["archives"]],
        )