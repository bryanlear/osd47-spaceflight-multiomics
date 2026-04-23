from pathlib import Path


output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

header = [
    "archive_id",
    "condition",
    "replicate",
    "lane_fastq",
    "lane_name",
]
lines = ["\t".join(header)]

extract_root = Path("results/extracted")

for entry in snakemake.config["archives"]:
    archive_dir = extract_root / entry["id"]
    fastqs = sorted(archive_dir.rglob("*.fastq.gz"))
    for fastq in fastqs:
        lines.append(
            "\t".join(
                [
                    entry["id"],
                    entry["condition"],
                    str(entry["replicate"]),
                    str(fastq),
                    fastq.name,
                ]
            )
        )

output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")