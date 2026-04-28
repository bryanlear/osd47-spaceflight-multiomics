from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests


EXPECTED_CONDITIONS = ("BSL", "FLT", "GC")
CONTRASTS = {
    "FLT_vs_BSL": {
        "left": "FLT",
        "right": "BSL",
        "hypothesis": "C(condition, Treatment(reference='BSL'))[T.FLT] = 0",
    },
    "GC_vs_BSL": {
        "left": "GC",
        "right": "BSL",
        "hypothesis": "C(condition, Treatment(reference='BSL'))[T.GC] = 0",
    },
    "FLT_vs_GC": {
        "left": "FLT",
        "right": "GC",
        "hypothesis": (
            "C(condition, Treatment(reference='BSL'))[T.FLT] - "
            "C(condition, Treatment(reference='BSL'))[T.GC] = 0"
        ),
    },
}


def default_output_dir() -> Path:
    return Path(__file__).resolve().parent / "results" / "fragpipe_tmt10_ms3_run3_protein_de"


def default_annotations_path() -> Path:
    return default_output_dir() / "protein_annotations.tsv"


def default_abundance_matrix_path() -> Path:
    return default_output_dir() / "protein_abundance_matrix.tsv"


def default_sample_metadata_path() -> Path:
    return default_output_dir() / "sample_metadata.tsv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Fit per-protein linear models for TMT protein abundance and test FLT-BSL, "
            "GC-BSL, and FLT-GC contrasts with BH FDR correction."
        )
    )
    parser.add_argument(
        "--annotations",
        type=Path,
        default=default_annotations_path(),
        help="Path to protein_annotations.tsv created by prepare_protein_de_inputs.py.",
    )
    parser.add_argument(
        "--abundance-matrix",
        type=Path,
        default=default_abundance_matrix_path(),
        help="Path to protein_abundance_matrix.tsv created by prepare_protein_de_inputs.py.",
    )
    parser.add_argument(
        "--sample-metadata",
        type=Path,
        default=default_sample_metadata_path(),
        help="Path to sample_metadata.tsv created by prepare_protein_de_inputs.py.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=default_output_dir(),
        help="Directory where result TSV files will be written.",
    )
    parser.add_argument(
        "--id-column",
        default="Protein",
        help="Protein identifier column shared across the annotation and abundance files.",
    )
    parser.add_argument(
        "--min-number-psm",
        type=float,
        default=2.0,
        help="Minimum NumberPSM required for a protein to be modeled.",
    )
    parser.add_argument(
        "--min-max-pep-prob",
        type=float,
        default=0.99,
        help="Minimum MaxPepProb required for a protein to be modeled.",
    )
    parser.add_argument(
        "--min-combined-spectral-count",
        type=float,
        default=0.0,
        help="Optional minimum Combined Spectral Count threshold.",
    )
    parser.add_argument(
        "--min-reps-per-condition",
        type=int,
        default=2,
        help="Minimum observed replicates required in each condition after NA removal.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="FDR threshold used to label significant proteins.",
    )
    return parser.parse_args()


def scalar(value: object) -> float:
    array = np.asarray(value)
    return float(array.reshape(-1)[0])


def bh_adjust(p_values: pd.Series) -> pd.Series:
    adjusted = pd.Series(np.nan, index=p_values.index, dtype=float)
    valid_mask = p_values.notna()
    if valid_mask.any():
        adjusted.loc[valid_mask] = multipletests(p_values.loc[valid_mask], method="fdr_bh")[1]
    return adjusted


def load_inputs(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str]]:
    annotations = pd.read_csv(args.annotations, sep="\t")
    abundance_matrix = pd.read_csv(args.abundance_matrix, sep="\t")
    sample_metadata = pd.read_csv(args.sample_metadata, sep="\t")

    for required in [args.id_column, "NumberPSM", "MaxPepProb"]:
        if required not in annotations.columns:
            raise ValueError(f"Required annotation column '{required}' was not found in {args.annotations}.")

    if args.id_column not in abundance_matrix.columns:
        raise ValueError(
            f"Identifier column '{args.id_column}' was not found in {args.abundance_matrix}."
        )

    sample_columns = sample_metadata["sample"].tolist()
    missing_samples = [sample for sample in sample_columns if sample not in abundance_matrix.columns]
    if missing_samples:
        raise ValueError(
            "Samples were present in sample_metadata.tsv but missing from the abundance matrix: "
            + ", ".join(missing_samples)
        )

    abundance_matrix = abundance_matrix.loc[:, [args.id_column, *sample_columns]].copy()
    abundance_matrix[sample_columns] = abundance_matrix[sample_columns].apply(
        pd.to_numeric,
        errors="coerce",
    )

    sample_metadata = sample_metadata.copy()
    sample_metadata["condition"] = pd.Categorical(
        sample_metadata["condition"],
        categories=list(EXPECTED_CONDITIONS),
        ordered=True,
    )

    return annotations, abundance_matrix, sample_metadata, sample_columns


def filter_annotations(annotations: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    filtered = annotations.copy()
    filtered["NumberPSM"] = pd.to_numeric(filtered["NumberPSM"], errors="coerce")
    filtered["MaxPepProb"] = pd.to_numeric(filtered["MaxPepProb"], errors="coerce")

    if "Combined Spectral Count" in filtered.columns:
        filtered["Combined Spectral Count"] = pd.to_numeric(
            filtered["Combined Spectral Count"],
            errors="coerce",
        )

    keep_mask = filtered["NumberPSM"].ge(args.min_number_psm) & filtered["MaxPepProb"].ge(
        args.min_max_pep_prob
    )
    if args.min_combined_spectral_count > 0:
        if "Combined Spectral Count" not in filtered.columns:
            raise ValueError(
                "--min-combined-spectral-count was set but Combined Spectral Count is absent "
                f"from {args.annotations}."
            )
        keep_mask &= filtered["Combined Spectral Count"].fillna(0).ge(args.min_combined_spectral_count)

    return filtered.loc[keep_mask].copy()


def build_long_table(
    abundance_matrix: pd.DataFrame,
    sample_metadata: pd.DataFrame,
    id_column: str,
    sample_columns: list[str],
) -> pd.DataFrame:
    long_table = abundance_matrix.melt(
        id_vars=[id_column],
        value_vars=sample_columns,
        var_name="sample",
        value_name="abundance",
    )
    long_table = long_table.merge(sample_metadata, on="sample", how="left", validate="many_to_one")
    long_table["abundance"] = pd.to_numeric(long_table["abundance"], errors="coerce")
    return long_table


def collect_results(
    filtered_annotations: pd.DataFrame,
    long_table: pd.DataFrame,
    args: argparse.Namespace,
    total_input_proteins: int,
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    annotation_lookup = filtered_annotations.set_index(args.id_column, drop=False)
    filtered_ids = annotation_lookup.index.tolist()
    filtered_long = long_table.loc[long_table[args.id_column].isin(filtered_ids)].dropna(
        subset=["abundance"]
    )

    contrast_rows: dict[str, list[dict[str, object]]] = {contrast: [] for contrast in CONTRASTS}
    modeled_ids: set[str] = set()
    skipped_low_rep_ids: set[str] = set()

    for protein_id, protein_table in filtered_long.groupby(args.id_column, sort=False):
        per_condition_counts = protein_table.groupby("condition", observed=False)["abundance"].count()
        if any(per_condition_counts.get(condition, 0) < args.min_reps_per_condition for condition in EXPECTED_CONDITIONS):
            skipped_low_rep_ids.add(str(protein_id))
            continue

        model = smf.ols(
            "abundance ~ C(condition, Treatment(reference='BSL'))",
            data=protein_table,
        ).fit()
        group_means = protein_table.groupby("condition", observed=False)["abundance"].mean()
        annotation_record = annotation_lookup.loc[protein_id].to_dict()
        modeled_ids.add(str(protein_id))

        for contrast_name, contrast_spec in CONTRASTS.items():
            test_result = model.t_test(contrast_spec["hypothesis"])
            left_condition = contrast_spec["left"]
            right_condition = contrast_spec["right"]

            record = dict(annotation_record)
            record.update(
                {
                    "contrast": contrast_name,
                    "left_condition": left_condition,
                    "right_condition": right_condition,
                    "estimate": scalar(test_result.effect),
                    "std_error": scalar(test_result.sd),
                    "t_stat": scalar(test_result.tvalue),
                    "p_value": scalar(test_result.pvalue),
                    "df_resid": float(model.df_resid),
                    "n_BSL": int(per_condition_counts.get("BSL", 0)),
                    "n_FLT": int(per_condition_counts.get("FLT", 0)),
                    "n_GC": int(per_condition_counts.get("GC", 0)),
                    "mean_BSL": float(group_means.get("BSL", np.nan)),
                    "mean_FLT": float(group_means.get("FLT", np.nan)),
                    "mean_GC": float(group_means.get("GC", np.nan)),
                    "mean_difference": float(
                        group_means.get(left_condition, np.nan) - group_means.get(right_condition, np.nan)
                    ),
                }
            )
            contrast_rows[contrast_name].append(record)

    contrast_tables: dict[str, pd.DataFrame] = {}
    summary_rows = [
        {
            "metric": "input_proteins",
            "value": int(total_input_proteins),
        },
        {
            "metric": "proteins_passing_confidence_filter",
            "value": int(len(filtered_annotations)),
        },
        {
            "metric": "proteins_filtered_out_confidence",
            "value": int(total_input_proteins - len(filtered_annotations)),
        },
        {
            "metric": "proteins_modeled",
            "value": int(len(modeled_ids)),
        },
        {
            "metric": "proteins_skipped_low_replicates",
            "value": int(len(skipped_low_rep_ids)),
        },
    ]

    for contrast_name, rows in contrast_rows.items():
        contrast_table = pd.DataFrame(rows)
        if contrast_table.empty:
            contrast_table["fdr_bh"] = pd.Series(dtype=float)
            contrast_table["significant"] = pd.Series(dtype=bool)
        else:
            contrast_table["fdr_bh"] = bh_adjust(contrast_table["p_value"])
            contrast_table["significant"] = contrast_table["fdr_bh"].le(args.alpha)
            contrast_table = contrast_table.sort_values(
                ["fdr_bh", "p_value", "estimate"],
                ascending=[True, True, False],
                na_position="last",
            )
        contrast_tables[contrast_name] = contrast_table
        summary_rows.append(
            {
                "metric": f"significant_{contrast_name}",
                "value": int(contrast_table["significant"].fillna(False).sum()) if not contrast_table.empty else 0,
            }
        )

    summary = pd.DataFrame(summary_rows)
    return contrast_tables, summary


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    annotations, abundance_matrix, sample_metadata, sample_columns = load_inputs(args)
    filtered_annotations = filter_annotations(annotations, args)
    long_table = build_long_table(abundance_matrix, sample_metadata, args.id_column, sample_columns)
    contrast_tables, summary = collect_results(
        filtered_annotations,
        long_table,
        args,
        total_input_proteins=len(annotations),
    )

    all_contrasts = []
    for contrast_name, contrast_table in contrast_tables.items():
        contrast_path = args.output_dir / f"differential_protein_abundance.{contrast_name}.tsv"
        contrast_table.to_csv(contrast_path, sep="\t", index=False)
        print(f"Wrote {contrast_name}: {contrast_path}")
        if not contrast_table.empty:
            all_contrasts.append(contrast_table)

    all_contrasts_path = args.output_dir / "differential_protein_abundance.all_contrasts.tsv"
    if all_contrasts:
        pd.concat(all_contrasts, axis=0, ignore_index=True).to_csv(
            all_contrasts_path,
            sep="\t",
            index=False,
        )
    else:
        pd.DataFrame().to_csv(all_contrasts_path, sep="\t", index=False)
    print(f"Wrote combined results: {all_contrasts_path}")

    summary_path = args.output_dir / "differential_protein_abundance.summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    print(f"Wrote summary: {summary_path}")


if __name__ == "__main__":
    main()