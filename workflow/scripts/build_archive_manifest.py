from pathlib import Path


output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

header = ["archive_id", "archive_file", "condition", "replicate", "exists"]
lines = ["\t".join(header)]

for entry in snakemake.config["archives"]:
    archive_file = Path(snakemake.config["input_dir"]) / entry["file"]
    lines.append(
        "\t".join(
            [
                entry["id"],
                str(archive_file),
                entry["condition"],
                str(entry["replicate"]),
                str(archive_file.exists()).lower(),
            ]
        )
    )

output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")