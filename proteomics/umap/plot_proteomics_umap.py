#!/usr/bin/env python3

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.preprocessing import StandardScaler
import umap


REPLICATE_PATTERN = re.compile(r"^(BSL|FLT|GC)_Rep(\d+)$")
COLORS = {
    "BSL": "#00c2a8",
    "FLT": "#ff6b35",
    "GC": "#3a86ff",
}
PANEL_BACKGROUND = "#ffffff"
UMAP_PARAMETERS_TEXT = "UMAP: metric=euclidean, n_neighbors=30, min_dist=0.15"


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def input_tables() -> dict[str, Path]:
    report_dir = repo_root() / "proteomics" / "results" / "fragpipe_tmt10_ms3_run3" / "tmt-report"
    return {
        "protein": report_dir / "abundance_protein_MD.tsv",
        "peptide": report_dir / "abundance_peptide_MD.tsv",
        "gene": report_dir / "abundance_gene_MD.tsv",
    }


def detect_condition_columns(table_path: Path) -> dict[str, list[str]]:
    header = pd.read_csv(table_path, sep="\t", nrows=0)
    condition_columns: dict[str, list[tuple[int, str]]] = {}

    for column in header.columns:
        match = REPLICATE_PATTERN.match(column)
        if not match:
            continue

        condition = match.group(1)
        replicate_number = int(match.group(2))
        condition_columns.setdefault(condition, []).append((replicate_number, column))

    ordered_columns = {
        condition: [column for _, column in sorted(replicates)]
        for condition, replicates in condition_columns.items()
    }

    if len(ordered_columns) < 2:
        raise ValueError(f"Could not find replicate columns in {table_path}")
    return ordered_columns


def load_condition_feature_matrix(table_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    header = pd.read_csv(table_path, sep="\t", nrows=0)
    id_column = "Index" if "Index" in header.columns else header.columns[0]
    condition_columns = detect_condition_columns(table_path)

    usecols = [id_column] + [column for columns in condition_columns.values() for column in columns]
    abundance = pd.read_csv(table_path, sep="\t", usecols=usecols)

    condition_frames = []
    for condition in ["BSL", "FLT", "GC"]:
        columns = condition_columns.get(condition)
        if not columns:
            continue

        condition_frame = abundance[[id_column] + columns].copy()
        replicate_columns = [f"Rep{index}" for index in range(1, len(columns) + 1)]
        condition_frame.columns = ["feature_id"] + replicate_columns
        condition_frame.insert(1, "condition", condition)
        condition_frames.append(condition_frame)

    feature_matrix = pd.concat(condition_frames, ignore_index=True)
    replicate_columns = [column for column in feature_matrix.columns if column.startswith("Rep")]

    feature_matrix[replicate_columns] = feature_matrix[replicate_columns].apply(pd.to_numeric, errors="coerce")
    feature_matrix = feature_matrix.dropna(subset=replicate_columns, how="all")

    row_medians = feature_matrix[replicate_columns].median(axis=1, skipna=True)
    feature_matrix[replicate_columns] = feature_matrix[replicate_columns].T.fillna(row_medians).T
    feature_matrix = feature_matrix.dropna(subset=replicate_columns, how="any")

    scaled_values = StandardScaler().fit_transform(feature_matrix[replicate_columns])
    scaled_matrix = pd.DataFrame(scaled_values, columns=replicate_columns, index=feature_matrix.index)
    metadata = feature_matrix[["feature_id", "condition"]].copy()

    return metadata, scaled_matrix


def build_embedding(feature_matrix: pd.DataFrame, random_state: int = 42) -> pd.DataFrame:
    n_samples = feature_matrix.shape[0]
    n_neighbors = min(30, n_samples - 1)
    if n_neighbors < 2:
        raise ValueError("UMAP requires at least 3 samples")

    reducer = umap.UMAP(
        n_components=2,
        n_neighbors=n_neighbors,
        min_dist=0.15,
        metric="euclidean",
        random_state=random_state,
        low_memory=True,
    )
    embedding = reducer.fit_transform(feature_matrix)

    return pd.DataFrame(embedding, columns=["UMAP1", "UMAP2"], index=feature_matrix.index)


def style_axis(ax: plt.Axes, x_limits: tuple[float, float], y_limits: tuple[float, float]) -> None:
    ax.set_facecolor(PANEL_BACKGROUND)
    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)
    ax.grid(False)
    ax.tick_params(labelsize=10)
    for spine in ax.spines.values():
        spine.set_color("#222222")
        spine.set_linewidth(0.8)


def plot_embedding(embedding: pd.DataFrame, title: str, output_prefix: Path) -> None:
    x_min, x_max = embedding["UMAP1"].min(), embedding["UMAP1"].max()
    y_min, y_max = embedding["UMAP2"].min(), embedding["UMAP2"].max()
    x_pad = (x_max - x_min) * 0.06
    y_pad = (y_max - y_min) * 0.06
    x_limits = (x_min - x_pad, x_max + x_pad)
    y_limits = (y_min - y_pad, y_max + y_pad)

    fig, ax = plt.subplots(figsize=(9.5, 7.5), dpi=150)
    fig.patch.set_facecolor(PANEL_BACKGROUND)

    style_axis(ax, x_limits, y_limits)

    for condition in ["BSL", "FLT", "GC"]:
        subset = embedding.loc[embedding["condition"] == condition]
        ax.scatter(
            subset["UMAP1"],
            subset["UMAP2"],
            color=COLORS[condition],
            label=f"{condition} (n={len(subset):,})",
            s=3.2,
            linewidth=0,
            alpha=0.72,
            rasterized=True,
        )
    ax.set_title(title, fontsize=16, fontweight="bold", pad=12)
    ax.set_xlabel("UMAP 1", fontsize=12)
    ax.set_ylabel("UMAP 2", fontsize=12)
    ax.legend(
        title="Condition",
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.0,
        frameon=True,
        facecolor="white",
        edgecolor="#d0d0d0",
        framealpha=0.95,
        markerscale=3.5,
    )
    fig.text(
        0.99,
        0.015,
        UMAP_PARAMETERS_TEXT,
        ha="right",
        va="bottom",
        fontsize=9,
        color="#6b7280",
    )
    fig.tight_layout(rect=(0.02, 0.04, 0.84, 0.97))
    fig.savefig(output_prefix.with_suffix(".png"))
    fig.savefig(output_prefix.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    output_dir = Path(__file__).resolve().parent
    output_dir.mkdir(parents=True, exist_ok=True)

    for label, table_path in input_tables().items():
        metadata, feature_matrix = load_condition_feature_matrix(table_path)
        embedding_coordinates = build_embedding(feature_matrix)
        embedding = pd.concat(
            [
                metadata.reset_index(drop=True),
                embedding_coordinates.reset_index(drop=True),
            ],
            axis=1,
        )

        output_prefix = output_dir / f"{label}_umap"
        plot_embedding(
            embedding,
            title=f"OSD-47 proteomics {label}-level UMAP by condition",
            output_prefix=output_prefix,
        )
        embedding.to_csv(output_prefix.with_suffix(".tsv"), sep="\t", index=False)

        counts = embedding["condition"].value_counts().reindex(["BSL", "FLT", "GC"])
        counts_summary = ", ".join(
            f"{condition}={int(count):,}"
            for condition, count in counts.items()
            if pd.notna(count)
        )

        print(
            f"{label}: {len(embedding):,} condition-specific points ({counts_summary}) -> "
            f"{output_prefix.with_suffix('.png').name}"
        )


if __name__ == "__main__":
    main()