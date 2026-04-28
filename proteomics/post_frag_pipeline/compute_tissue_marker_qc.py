from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


MUSCLE_MARKERS = [
    "Myh1",
    "Myh2",
    "Myh4",
    "Myh7",
    "Myh7b",
    "Acta1",
    "Tnnt1",
    "Tnnt3",
    "Tnni2",
    "Myl1",
    "Myl2",
    "Myl3",
    "Ckm",
    "Mb",
    "Casq1",
    "Pvalb",
]

LIVER_MARKERS = [
    "Alb",
    "Ttr",
    "Apoa1",
    "Apoa2",
    "Apob",
    "Cps1",
    "Ass1",
    "Asl",
    "Arg1",
    "Cyp2e1",
    "Cyp3a11",
    "Fasn",
    "Acaca",
    "Pck1",
    "G6pc",
]

GENE_SYMBOL_PATTERN = re.compile(r"gene_symbol:([^ ]+)")


def default_results_dir() -> Path:
    return Path(__file__).resolve().parent / "results" / "fragpipe_tmt10_ms3_run3_protein_de"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute per-sample tissue-marker QC scores from the prepared post-FragPipe protein abundance tables."
        )
    )
    parser.add_argument(
        "--annotations",
        type=Path,
        default=default_results_dir() / "protein_annotations.tsv",
        help="Prepared protein annotation table.",
    )
    parser.add_argument(
        "--abundance-matrix",
        type=Path,
        default=default_results_dir() / "protein_abundance_matrix.tsv",
        help="Prepared protein abundance matrix containing only biological replicates.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=default_results_dir(),
        help="Directory where tissue marker QC outputs will be written.",
    )
    parser.add_argument(
        "--id-column",
        default="Protein",
        help="Protein identifier column shared between annotation and abundance tables.",
    )
    return parser.parse_args()


def parse_gene_symbol(row: pd.Series) -> str:
    gene = row.get("Gene")
    if isinstance(gene, str) and gene.strip():
        return gene.strip()

    for column in ["Protein", "Protein ID", "Entry Name", "Protein Description", "Index"]:
        value = row.get(column)
        if not isinstance(value, str):
            continue
        match = GENE_SYMBOL_PATTERN.search(value)
        if match is not None:
            return match.group(1)
    return ""


def classify_panel(gene_symbol: str) -> list[str]:
    panels: list[str] = []
    if gene_symbol in MUSCLE_MARKERS:
        panels.append("muscle")
    if gene_symbol in LIVER_MARKERS or gene_symbol.startswith("Cyp2c"):
        panels.append("liver")
    return panels


def zscore_series(values: pd.Series) -> pd.Series:
    std = values.std(ddof=1)
    if pd.isna(std) or std == 0:
        return pd.Series(0.0, index=values.index)
    return (values - values.mean()) / std


def choose_marker_representatives(
    annotations: pd.DataFrame,
    abundance_matrix: pd.DataFrame,
    id_column: str,
) -> tuple[pd.DataFrame, list[str]]:
    sample_columns = [column for column in abundance_matrix.columns if column != id_column]

    annotations = annotations.copy()
    abundance_matrix = abundance_matrix.copy()
    annotations[id_column] = annotations[id_column].astype(str)
    abundance_matrix[id_column] = abundance_matrix[id_column].astype(str)

    merged = annotations.merge(abundance_matrix, on=id_column, how="inner", validate="one_to_one")
    merged["gene_symbol"] = merged.apply(parse_gene_symbol, axis=1)
    merged["panel_list"] = merged["gene_symbol"].map(classify_panel)
    merged = merged.loc[merged["panel_list"].map(bool)].copy()

    if merged.empty:
        raise ValueError("No requested tissue markers were found in the prepared protein tables.")

    merged["Combined Spectral Count"] = pd.to_numeric(merged.get("Combined Spectral Count"), errors="coerce")
    merged["NumberPSM"] = pd.to_numeric(merged.get("NumberPSM"), errors="coerce")
    merged["MaxPepProb"] = pd.to_numeric(merged.get("MaxPepProb"), errors="coerce")
    merged[sample_columns] = merged[sample_columns].apply(pd.to_numeric, errors="coerce")

    representatives: list[pd.Series] = []
    for gene_symbol, subset in merged.groupby("gene_symbol", sort=True):
        ranked = subset.sort_values(
            ["Combined Spectral Count", "NumberPSM", "MaxPepProb"],
            ascending=[False, False, False],
            na_position="last",
        )
        representatives.append(ranked.iloc[0])

    representative_df = pd.DataFrame(representatives).reset_index(drop=True)
    representative_df["panel"] = representative_df["panel_list"].map(lambda panels: ",".join(panels))
    return representative_df, sample_columns


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    annotations = pd.read_csv(args.annotations, sep="\t")
    abundance_matrix = pd.read_csv(args.abundance_matrix, sep="\t")

    representative_df, sample_columns = choose_marker_representatives(
        annotations,
        abundance_matrix,
        args.id_column,
    )

    marker_selection = representative_df[
        [
            "panel",
            "gene_symbol",
            args.id_column,
            "Combined Spectral Count",
            "NumberPSM",
            "MaxPepProb",
            *sample_columns,
        ]
    ].copy()
    marker_selection = marker_selection.sort_values(["panel", "gene_symbol"]).reset_index(drop=True)

    zscore_table = marker_selection[["panel", "gene_symbol", args.id_column]].copy()
    zscore_values = marker_selection.loc[:, sample_columns].apply(zscore_series, axis=1)
    zscore_values.index = marker_selection.index
    zscore_table = pd.concat([zscore_table, zscore_values], axis=1)

    score_rows = []
    for panel_name in ["muscle", "liver"]:
        panel_mask = zscore_table["panel"].str.contains(panel_name)
        panel_table = zscore_table.loc[panel_mask, sample_columns]
        if panel_table.empty:
            continue
        panel_mean = panel_table.mean(axis=0)
        for sample, score in panel_mean.items():
            score_rows.append(
                {
                    "sample": sample,
                    "panel": panel_name,
                    "mean_zscore": float(score),
                    "n_markers_used": int(panel_table.shape[0]),
                }
            )

    sample_scores = pd.DataFrame(score_rows)
    wide_scores = sample_scores.pivot(index="sample", columns="panel", values="mean_zscore").reset_index()
    wide_scores["muscle_minus_liver"] = wide_scores.get("muscle", 0.0) - wide_scores.get("liver", 0.0)
    wide_scores = wide_scores.sort_values("sample").reset_index(drop=True)

    marker_selection_path = args.output_dir / "tissue_marker_selection.tsv"
    gene_zscores_path = args.output_dir / "tissue_marker_gene_zscores.tsv"
    sample_scores_path = args.output_dir / "tissue_marker_sample_scores.tsv"

    marker_selection.to_csv(marker_selection_path, sep="\t", index=False)
    zscore_table.to_csv(gene_zscores_path, sep="\t", index=False)
    wide_scores.to_csv(sample_scores_path, sep="\t", index=False)

    print(f"Wrote marker selection: {marker_selection_path}")
    print(f"Wrote marker gene z-scores: {gene_zscores_path}")
    print(f"Wrote sample tissue scores: {sample_scores_path}")
    print(wide_scores.to_csv(sep='\t', index=False))


if __name__ == "__main__":
    main()