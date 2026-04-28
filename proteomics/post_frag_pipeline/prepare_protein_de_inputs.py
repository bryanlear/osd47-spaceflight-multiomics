from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


EXPECTED_CONDITIONS = ("BSL", "FLT", "GC")
REPLICATE_PATTERN = re.compile(r"^(BSL|FLT|GC)_Rep(\d+)$")


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def default_abundance_path() -> Path:
    return (
        repo_root()
        / "proteomics"
        / "results"
        / "fragpipe_tmt10_ms3_run3"
        / "tmt-report"
        / "abundance_protein_MD.tsv"
    )


def default_combined_path() -> Path:
    return (
        repo_root()
        / "proteomics"
        / "results"
        / "fragpipe_tmt10_ms3_run3"
        / "combined_protein.tsv"
    )


def default_output_dir() -> Path:
    return Path(__file__).resolve().parent / "results" / "fragpipe_tmt10_ms3_run3_protein_de"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Split protein-level annotations from the 9 biological replicate abundance columns "
            "and join Combined Spectral Count from the FragPipe combined_protein table."
        )
    )
    parser.add_argument(
        "--abundance",
        type=Path,
        default=default_abundance_path(),
        help="Path to abundance_protein_MD.tsv.",
    )
    parser.add_argument(
        "--combined",
        type=Path,
        default=default_combined_path(),
        help="Path to combined_protein.tsv.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=default_output_dir(),
        help="Directory where the prepared TSV files will be written.",
    )
    parser.add_argument(
        "--join-key",
        default="Protein",
        choices=["Protein", "Protein ID"],
        help="Column used to join abundance_protein_MD.tsv to combined_protein.tsv.",
    )
    return parser.parse_args()


def ordered_replicate_columns(columns: list[str]) -> list[str]:
    matched_columns: list[tuple[int, int, str]] = []
    condition_order = {condition: index for index, condition in enumerate(EXPECTED_CONDITIONS)}

    for column in columns:
        match = REPLICATE_PATTERN.match(column)
        if match is None:
            continue
        condition, replicate = match.groups()
        matched_columns.append((condition_order[condition], int(replicate), column))

    if not matched_columns:
        raise ValueError("No biological replicate columns matching BSL/FLT/GC were found.")

    matched_columns.sort()
    replicate_columns = [column for _, _, column in matched_columns]

    observed_conditions = {column.split("_Rep", maxsplit=1)[0] for column in replicate_columns}
    missing_conditions = [condition for condition in EXPECTED_CONDITIONS if condition not in observed_conditions]
    if missing_conditions:
        raise ValueError(
            "Missing replicate columns for conditions: " + ", ".join(missing_conditions)
        )

    return replicate_columns


def load_combined_spectral_counts(combined_path: Path, join_key: str) -> pd.DataFrame:
    combined = pd.read_csv(combined_path, sep="\t", usecols=[join_key, "Combined Spectral Count"])
    combined = combined.dropna(subset=[join_key]).drop_duplicates(subset=[join_key], keep="first")
    return combined


def build_sample_metadata(replicate_columns: list[str]) -> pd.DataFrame:
    rows = []
    for sample in replicate_columns:
        condition, replicate = sample.split("_Rep", maxsplit=1)
        rows.append(
            {
                "sample": sample,
                "condition": condition,
                "replicate": int(replicate),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    abundance = pd.read_csv(args.abundance, sep="\t")
    if args.join_key not in abundance.columns:
        raise ValueError(f"Join key '{args.join_key}' was not found in {args.abundance}.")

    if abundance[args.join_key].duplicated().any():
        raise ValueError(f"Join key '{args.join_key}' is not unique in {args.abundance}.")

    replicate_columns = ordered_replicate_columns(abundance.columns.tolist())
    combined_spectral_counts = load_combined_spectral_counts(args.combined, args.join_key)

    prepared = abundance.merge(
        combined_spectral_counts,
        on=args.join_key,
        how="left",
        validate="one_to_one",
    )

    annotation_columns = [column for column in prepared.columns if column not in replicate_columns]
    annotations = prepared.loc[:, annotation_columns].copy()
    abundance_matrix = prepared.loc[:, [args.join_key, *replicate_columns]].copy()
    sample_metadata = build_sample_metadata(replicate_columns)

    annotations_path = args.output_dir / "protein_annotations.tsv"
    abundance_matrix_path = args.output_dir / "protein_abundance_matrix.tsv"
    sample_metadata_path = args.output_dir / "sample_metadata.tsv"

    annotations.to_csv(annotations_path, sep="\t", index=False)
    abundance_matrix.to_csv(abundance_matrix_path, sep="\t", index=False)
    sample_metadata.to_csv(sample_metadata_path, sep="\t", index=False)

    unmatched = int(prepared["Combined Spectral Count"].isna().sum())
    print(f"Wrote annotations: {annotations_path}")
    print(f"Wrote abundance matrix: {abundance_matrix_path}")
    print(f"Wrote sample metadata: {sample_metadata_path}")
    print(f"Prepared proteins: {len(prepared)}")
    print(f"Replicate columns: {', '.join(replicate_columns)}")
    print(f"Proteins without Combined Spectral Count: {unmatched}")


if __name__ == "__main__":
    main()