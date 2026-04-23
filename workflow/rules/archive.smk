from pathlib import Path


ARCHIVE_BY_ID = {entry["id"]: entry for entry in config["archives"]}


def archive_path(archive_id):
    entry = ARCHIVE_BY_ID[archive_id]
    return str(Path(config["input_dir"]) / entry["file"])


rule build_archive_manifest:
    input:
        lambda wildcards: [archive_path(entry["id"]) for entry in config["archives"]]
    output:
        config["outputs"]["archive_manifest"]
    script:
        "../scripts/build_archive_manifest.py"


rule extract_archive:
    input:
        archive=lambda wildcards: archive_path(wildcards.archive)
    output:
        marker="results/extracted/{archive}/.extracted.ok"
    params:
        outdir=lambda wildcards: f"results/extracted/{wildcards.archive}"
    shell:
        """
        mkdir -p {params.outdir}
        tar -xzf {input.archive} -C {params.outdir}
        touch {output.marker}
        """


rule build_lane_manifest:
    input:
        markers=expand(
            "results/extracted/{archive}/.extracted.ok",
            archive=[entry["id"] for entry in config["archives"]],
        )
    output:
        config["outputs"]["lane_manifest"]
    script:
        "../scripts/build_lane_manifest.py"