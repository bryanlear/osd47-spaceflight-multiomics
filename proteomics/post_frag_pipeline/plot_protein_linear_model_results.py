from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D


CONDITION_PALETTE = {
    "BSL": "#0072B2",
    "FLT": "#D55E00",
    "GC": "#009E73",
}

PLOT_NEUTRALS = {
    "ink": "#222222",
    "grid": "#D9E2EC",
    "nonsig": "#C7CDD4",
    "sig_up": "#C45A1A",
    "sig_down": "#3A6EA5",
    "threshold": "#6B7280",
    "threshold_light": "#9CA3AF",
}

PCA_LABEL_OFFSETS = [(4, 4), (4, -10), (6, 10), (-28, 4), (6, -12)]
VOLCANO_POSITIVE_LABEL_OFFSETS = [(8, 8), (8, -12), (10, 14), (12, -18), (6, 20)]
VOLCANO_NEGATIVE_LABEL_OFFSETS = [(-8, 8), (-8, -12), (-10, 14), (-12, -18), (-6, 20)]
GENE_SYMBOL_PATTERN = re.compile(r"gene_symbol:([^ ]+)")
EXPECTED_CONDITION_ORDER = ["BSL", "FLT", "GC"]
CONTRAST_FILES = ["FLT_vs_BSL", "GC_vs_BSL", "FLT_vs_GC"]


def default_output_dir() -> Path:
    return Path(__file__).resolve().parent / "results" / "fragpipe_tmt10_ms3_run3_protein_de"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create protein-level PCA, heatmap, and volcano plots from the post-FragPipe linear-model outputs."
        )
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=default_output_dir(),
        help="Directory containing the protein DE result and intermediate TSV files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=default_output_dir(),
        help="Directory where the plots and plot-sidecar TSV files will be written.",
    )
    parser.add_argument(
        "--id-column",
        default="Protein",
        help="Protein identifier column shared across the annotations and abundance matrix.",
    )
    parser.add_argument(
        "--min-number-psm",
        type=float,
        default=2.0,
        help="Minimum NumberPSM required for proteins to contribute to PCA and heatmaps.",
    )
    parser.add_argument(
        "--min-max-pep-prob",
        type=float,
        default=0.99,
        help="Minimum MaxPepProb required for proteins to contribute to PCA and heatmaps.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="Adjusted p-value threshold used to color volcano plots.",
    )
    parser.add_argument(
        "--effect-threshold",
        type=float,
        default=0.5,
        help="Absolute abundance-difference threshold to annotate on volcano plots.",
    )
    parser.add_argument(
        "--top-labels-per-direction",
        type=int,
        default=5,
        help="Number of top positive and negative proteins to label per volcano plot.",
    )
    parser.add_argument(
        "--heatmap-proteins",
        type=int,
        default=50,
        help="Number of highest-variance proteins to include in the protein heatmap.",
    )
    return parser.parse_args()


def style_axis(ax: plt.Axes, grid_axis: str = "y") -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(PLOT_NEUTRALS["ink"])
    ax.spines["bottom"].set_color(PLOT_NEUTRALS["ink"])
    ax.tick_params(colors=PLOT_NEUTRALS["ink"])
    if grid_axis:
        ax.grid(True, axis=grid_axis, color=PLOT_NEUTRALS["grid"], linewidth=0.8)
    ax.set_axisbelow(True)


def save_figure(fig: plt.Figure, output_path: Path) -> None:
    fig.savefig(output_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def contrast_label(contrast_name: str) -> str:
    left, right = contrast_name.split("_vs_")
    return f"{left} vs {right}"


def parse_gene_symbol(value: object) -> str | None:
    if not isinstance(value, str) or not value.strip():
        return None
    match = GENE_SYMBOL_PATTERN.search(value)
    if match is None:
        return None
    return match.group(1)


def protein_label(row: pd.Series, id_column: str) -> str:
    gene_value = row.get("Gene")
    if isinstance(gene_value, str) and gene_value.strip():
        return gene_value.strip()

    for column in ["Protein", "Protein ID", "Entry Name", "Protein Description"]:
        parsed = parse_gene_symbol(row.get(column))
        if parsed:
            return parsed

    protein_id = row.get(id_column)
    if isinstance(protein_id, str) and protein_id.strip():
        return protein_id.split()[0]
    return str(protein_id)


def load_inputs(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    annotations = pd.read_csv(args.input_dir / "protein_annotations.tsv", sep="\t")
    abundance_matrix = pd.read_csv(args.input_dir / "protein_abundance_matrix.tsv", sep="\t")
    sample_metadata = pd.read_csv(args.input_dir / "sample_metadata.tsv", sep="\t")
    all_contrasts = pd.read_csv(args.input_dir / "differential_protein_abundance.all_contrasts.tsv", sep="\t")
    return annotations, abundance_matrix, sample_metadata, all_contrasts


def filter_plot_matrix(
    annotations: pd.DataFrame,
    abundance_matrix: pd.DataFrame,
    sample_metadata: pd.DataFrame,
    modeled_ids: set[str],
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    annotations = annotations.copy()
    annotations["NumberPSM"] = pd.to_numeric(annotations["NumberPSM"], errors="coerce")
    annotations["MaxPepProb"] = pd.to_numeric(annotations["MaxPepProb"], errors="coerce")

    keep_mask = annotations[args.id_column].isin(modeled_ids)
    keep_mask &= annotations["NumberPSM"].ge(args.min_number_psm)
    keep_mask &= annotations["MaxPepProb"].ge(args.min_max_pep_prob)
    filtered_annotations = annotations.loc[keep_mask].copy()

    sample_metadata = sample_metadata.copy()
    sample_metadata["condition"] = pd.Categorical(
        sample_metadata["condition"],
        categories=EXPECTED_CONDITION_ORDER,
        ordered=True,
    )
    sample_metadata = sample_metadata.sort_values(["condition", "replicate"]).reset_index(drop=True)
    sample_columns = sample_metadata["sample"].tolist()

    filtered_matrix = abundance_matrix.loc[
        abundance_matrix[args.id_column].isin(filtered_annotations[args.id_column]),
        [args.id_column, *sample_columns],
    ].copy()
    filtered_matrix = filtered_matrix.drop_duplicates(subset=[args.id_column], keep="first")
    filtered_matrix = filtered_matrix.set_index(args.id_column)
    filtered_matrix = filtered_matrix.apply(pd.to_numeric, errors="coerce")

    complete_case_matrix = filtered_matrix.dropna(axis=0, how="any")
    if len(complete_case_matrix.index) >= 2:
        return complete_case_matrix, sample_metadata, "complete_cases"

    imputed_matrix = filtered_matrix.apply(lambda row: row.fillna(row.mean()), axis=1)
    imputed_matrix = imputed_matrix.dropna(axis=0, how="any")
    if len(imputed_matrix.index) < 2:
        raise ValueError("Not enough proteins remain to build PCA and heatmap plots.")
    return imputed_matrix, sample_metadata, "row_mean_imputed"


def select_volcano_labels(result: pd.DataFrame, alpha: float, labels_per_direction: int) -> pd.DataFrame:
    candidates = result.loc[result["fdr_bh"].notna() & (result["fdr_bh"] < alpha)].copy()
    if candidates.empty:
        return candidates

    candidates["abs_estimate"] = candidates["estimate"].abs()
    up = candidates.loc[candidates["estimate"] > 0].sort_values(
        ["fdr_bh", "abs_estimate"],
        ascending=[True, False],
    ).head(labels_per_direction)
    down = candidates.loc[candidates["estimate"] < 0].sort_values(
        ["fdr_bh", "abs_estimate"],
        ascending=[True, False],
    ).head(labels_per_direction)
    return pd.concat([up, down], axis=0).sort_values(["fdr_bh", "abs_estimate"], ascending=[True, False])


def annotate_volcano_labels(ax: plt.Axes, labels: pd.DataFrame) -> None:
    positive_index = 0
    negative_index = 0
    for row in labels.itertuples(index=False):
        if row.estimate >= 0:
            dx, dy = VOLCANO_POSITIVE_LABEL_OFFSETS[positive_index % len(VOLCANO_POSITIVE_LABEL_OFFSETS)]
            positive_index += 1
        else:
            dx, dy = VOLCANO_NEGATIVE_LABEL_OFFSETS[negative_index % len(VOLCANO_NEGATIVE_LABEL_OFFSETS)]
            negative_index += 1

        ax.annotate(
            row.plot_label,
            (row.estimate, row.neg_log10_fdr),
            xytext=(dx, dy),
            textcoords="offset points",
            ha="left" if row.estimate >= 0 else "right",
            va="center",
            fontsize=7.5,
            color=PLOT_NEUTRALS["ink"],
            bbox={"boxstyle": "round,pad=0.18", "facecolor": "white", "edgecolor": "none", "alpha": 0.82},
            arrowprops={"arrowstyle": "-", "color": PLOT_NEUTRALS["threshold_light"], "linewidth": 0.7},
            clip_on=False,
        )


def save_volcano_plots(annotations: pd.DataFrame, args: argparse.Namespace) -> None:
    annotation_lookup = annotations.set_index(args.id_column, drop=False)

    for contrast_name in CONTRAST_FILES:
        result_path = args.input_dir / f"differential_protein_abundance.{contrast_name}.tsv"
        result = pd.read_csv(result_path, sep="\t").copy()

        if args.id_column not in result.columns:
            raise ValueError(f"Identifier column '{args.id_column}' is missing from {result_path}.")

        result["plot_label"] = result[args.id_column].map(
            lambda protein_id: protein_label(annotation_lookup.loc[protein_id], args.id_column)
            if protein_id in annotation_lookup.index
            else str(protein_id)
        )
        result["estimate"] = pd.to_numeric(result["estimate"], errors="coerce")
        result["fdr_bh"] = pd.to_numeric(result["fdr_bh"], errors="coerce")
        fdr_for_plot = result["fdr_bh"].replace(0, np.nextafter(0, 1))
        neg_log10_fdr = -np.log10(fdr_for_plot)

        significant_mask = result["fdr_bh"].notna() & (result["fdr_bh"] < args.alpha)
        up = result.loc[significant_mask & (result["estimate"] > 0)].copy()
        down = result.loc[significant_mask & (result["estimate"] < 0)].copy()
        nonsig = result.loc[~significant_mask].copy()

        labels = select_volcano_labels(result, args.alpha, args.top_labels_per_direction)
        if not labels.empty:
            labels = labels.copy()
            labels["neg_log10_fdr"] = -np.log10(labels["fdr_bh"].replace(0, np.nextafter(0, 1)))

        fig, ax = plt.subplots(figsize=(8.2, 5.8))
        ax.scatter(
            nonsig["estimate"],
            neg_log10_fdr.loc[nonsig.index],
            s=12,
            c=PLOT_NEUTRALS["nonsig"],
            alpha=0.65,
        )
        ax.scatter(
            up["estimate"],
            -np.log10(up["fdr_bh"].replace(0, np.nextafter(0, 1))),
            s=18,
            c=PLOT_NEUTRALS["sig_up"],
            alpha=0.82,
        )
        ax.scatter(
            down["estimate"],
            -np.log10(down["fdr_bh"].replace(0, np.nextafter(0, 1))),
            s=18,
            c=PLOT_NEUTRALS["sig_down"],
            alpha=0.82,
        )
        ax.axhline(-np.log10(args.alpha), color=PLOT_NEUTRALS["threshold"], linestyle="--", linewidth=1)
        ax.axvline(args.effect_threshold, color=PLOT_NEUTRALS["threshold_light"], linestyle=":", linewidth=1)
        ax.axvline(-args.effect_threshold, color=PLOT_NEUTRALS["threshold_light"], linestyle=":", linewidth=1)
        ax.set_xlabel("Estimated abundance difference")
        ax.set_ylabel("-log10 adjusted p-value")
        ax.set_title(f"Volcano plot: {contrast_label(contrast_name)}")
        style_axis(ax, grid_axis="y")
        if not labels.empty:
            annotate_volcano_labels(ax, labels)
        ax.legend(
            handles=[
                Line2D([0], [0], marker="o", linestyle="", color=PLOT_NEUTRALS["nonsig"], label="Not significant", markersize=6),
                Line2D([0], [0], marker="o", linestyle="", color=PLOT_NEUTRALS["sig_up"], label="Significant up", markersize=6),
                Line2D([0], [0], marker="o", linestyle="", color=PLOT_NEUTRALS["sig_down"], label="Significant down", markersize=6),
                Line2D([0], [0], linestyle="--", color=PLOT_NEUTRALS["threshold"], label=f"FDR = {args.alpha:g}"),
                Line2D([0], [0], linestyle=":", color=PLOT_NEUTRALS["threshold_light"], label=f"|effect| = {args.effect_threshold:g}"),
            ],
            frameon=False,
            loc="upper left",
            bbox_to_anchor=(1.01, 1.0),
            borderaxespad=0.0,
            title="Legend",
        )
        fig.subplots_adjust(right=0.78)
        save_figure(fig, args.output_dir / f"differential_protein_abundance.{contrast_name}.volcano.png")


def compute_pca(sample_matrix: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    centered = sample_matrix - sample_matrix.mean(axis=0)
    nonconstant = centered.loc[:, centered.var(axis=0, ddof=1).fillna(0) > 0]
    if nonconstant.shape[1] < 2:
        raise ValueError("Not enough variable proteins remain after filtering to compute PCA.")

    values = nonconstant.to_numpy(dtype=float)
    u, singular_values, _ = np.linalg.svd(values, full_matrices=False)
    coords = u[:, :2] * singular_values[:2]

    denominator = max(values.shape[0] - 1, 1)
    variance = (singular_values ** 2) / denominator
    explained_variance_ratio = variance / variance.sum()

    pca_df = pd.DataFrame({"PC1": coords[:, 0], "PC2": coords[:, 1]}, index=sample_matrix.index)
    return pca_df, explained_variance_ratio


def save_pca_plot(plot_matrix: pd.DataFrame, sample_metadata: pd.DataFrame, args: argparse.Namespace) -> None:
    sample_matrix = plot_matrix.transpose()
    pca_df, explained_variance_ratio = compute_pca(sample_matrix)
    metadata_lookup = sample_metadata.set_index("sample")
    pca_df["condition"] = pca_df.index.map(metadata_lookup["condition"])
    pca_df["replicate"] = pca_df.index.map(metadata_lookup["replicate"])
    pca_df.to_csv(args.output_dir / "pca_coordinates.tsv", sep="\t")

    fig, ax = plt.subplots(figsize=(8.2, 6.0))
    style_axis(ax, grid_axis="both")
    for condition, subset in pca_df.groupby("condition", observed=False):
        ax.scatter(
            subset["PC1"],
            subset["PC2"],
            label=condition,
            color=CONDITION_PALETTE[str(condition)],
            s=95,
            edgecolors="white",
            linewidth=0.8,
        )
        for index, (sample, row) in enumerate(subset.iterrows()):
            dx, dy = PCA_LABEL_OFFSETS[index % len(PCA_LABEL_OFFSETS)]
            ax.annotate(
                sample,
                (row["PC1"], row["PC2"]),
                xytext=(dx, dy),
                textcoords="offset points",
                fontsize=8,
                color=PLOT_NEUTRALS["ink"],
                clip_on=False,
            )

    ax.set_xlabel(f"PC1 ({explained_variance_ratio[0] * 100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({explained_variance_ratio[1] * 100:.1f}% variance)")
    ax.set_title("PCA on modeled protein abundances")
    ax.margins(x=0.08, y=0.08)
    ax.legend(title="Condition", frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)
    fig.subplots_adjust(right=0.82)
    save_figure(fig, args.output_dir / "pca_plot.png")


def save_sample_correlation_heatmap(plot_matrix: pd.DataFrame, sample_metadata: pd.DataFrame, args: argparse.Namespace) -> None:
    sample_matrix = plot_matrix.transpose()
    ordered_samples = sample_metadata["sample"].tolist()
    correlation = sample_matrix.transpose().corr().loc[ordered_samples, ordered_samples]
    correlation.to_csv(args.output_dir / "sample_correlation_matrix.tsv", sep="\t")

    cmap = LinearSegmentedColormap.from_list(
        "sample_corr",
        ["#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"],
    )
    fig, ax = plt.subplots(figsize=(8.5, 6.8))
    image = ax.imshow(correlation.values, cmap=cmap, vmin=correlation.values.min(), vmax=1.0)
    ax.set_xticks(range(len(ordered_samples)))
    ax.set_yticks(range(len(ordered_samples)))
    ax.set_xticklabels(ordered_samples, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_yticklabels(ordered_samples)
    ax.set_title("Sample-sample Pearson correlation")
    for idx, sample in enumerate(ordered_samples):
        condition = str(sample_metadata.loc[sample_metadata["sample"] == sample, "condition"].iloc[0])
        ax.get_xticklabels()[idx].set_color(CONDITION_PALETTE[condition])
        ax.get_yticklabels()[idx].set_color(CONDITION_PALETTE[condition])
    fig.colorbar(image, ax=ax, label="Pearson r", fraction=0.046, pad=0.04)
    save_figure(fig, args.output_dir / "sample_correlation_heatmap.png")


def save_top_variable_heatmap(
    plot_matrix: pd.DataFrame,
    annotations: pd.DataFrame,
    sample_metadata: pd.DataFrame,
    args: argparse.Namespace,
) -> None:
    variances = plot_matrix.var(axis=1, ddof=1).sort_values(ascending=False)
    top_ids = variances.head(args.heatmap_proteins).index.tolist()
    heatmap_matrix = plot_matrix.loc[top_ids, sample_metadata["sample"].tolist()].copy()
    row_means = heatmap_matrix.mean(axis=1)
    row_stds = heatmap_matrix.std(axis=1, ddof=1).replace(0, 1)
    zscore_matrix = heatmap_matrix.sub(row_means, axis=0).div(row_stds, axis=0)

    annotation_lookup = annotations.set_index(args.id_column, drop=False)
    display_labels = [protein_label(annotation_lookup.loc[protein_id], args.id_column) for protein_id in top_ids]

    zscore_matrix.to_csv(args.output_dir / "top_variable_proteins.zscore_matrix.tsv", sep="\t")

    fig, ax = plt.subplots(figsize=(8.8, max(7.0, len(top_ids) * 0.22)))
    image = ax.imshow(zscore_matrix.values, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    ax.set_xticks(range(len(zscore_matrix.columns)))
    ax.set_xticklabels(zscore_matrix.columns, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_yticks(range(len(display_labels)))
    ax.set_yticklabels(display_labels)
    ax.set_title(f"Top {len(top_ids)} most variable modeled proteins")
    for idx, sample in enumerate(zscore_matrix.columns):
        condition = str(sample_metadata.loc[sample_metadata["sample"] == sample, "condition"].iloc[0])
        ax.get_xticklabels()[idx].set_color(CONDITION_PALETTE[condition])
    fig.colorbar(image, ax=ax, label="Row z-score", fraction=0.046, pad=0.04)
    save_figure(fig, args.output_dir / "top_variable_proteins_heatmap.png")


def save_plot_summary(plot_matrix: pd.DataFrame, matrix_mode: str, sample_metadata: pd.DataFrame, args: argparse.Namespace) -> None:
    summary = pd.DataFrame(
        [
            {"metric": "plot_matrix_mode", "value": matrix_mode},
            {"metric": "proteins_in_plot_matrix", "value": int(len(plot_matrix.index))},
            {"metric": "samples_in_plot_matrix", "value": int(len(sample_metadata.index))},
            {"metric": "heatmap_proteins", "value": int(min(args.heatmap_proteins, len(plot_matrix.index)))},
        ]
    )
    summary.to_csv(args.output_dir / "plot_summary.tsv", sep="\t", index=False)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    annotations, abundance_matrix, sample_metadata, all_contrasts = load_inputs(args)
    modeled_ids = set(all_contrasts[args.id_column].dropna().astype(str))
    abundance_matrix[args.id_column] = abundance_matrix[args.id_column].astype(str)
    annotations[args.id_column] = annotations[args.id_column].astype(str)

    plot_matrix, ordered_metadata, matrix_mode = filter_plot_matrix(
        annotations,
        abundance_matrix,
        sample_metadata,
        modeled_ids,
        args,
    )
    plot_matrix.to_csv(args.output_dir / f"protein_plot_matrix.{matrix_mode}.tsv", sep="\t")

    save_volcano_plots(annotations, args)
    save_pca_plot(plot_matrix, ordered_metadata, args)
    save_sample_correlation_heatmap(plot_matrix, ordered_metadata, args)
    save_top_variable_heatmap(plot_matrix, annotations, ordered_metadata, args)
    save_plot_summary(plot_matrix, matrix_mode, ordered_metadata, args)

    print(f"Wrote volcano, PCA, and heatmap plots to {args.output_dir}")


if __name__ == "__main__":
    main()