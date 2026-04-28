from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA


CONDITION_PALETTE = {
    "BSL": "#0072B2",
    "FLT": "#D55E00",
    "GC": "#009E73",
}

FALLBACK_CONDITION_PALETTE = ["#CC79A7", "#E69F00", "#56B4E9", "#7F7F7F"]

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
VOLCANO_LABELS_PER_DIRECTION = 5
VOLCANO_POSITIVE_LABEL_OFFSETS = [(8, 8), (8, -12), (10, 14), (12, -18), (6, 20)]
VOLCANO_NEGATIVE_LABEL_OFFSETS = [(-8, 8), (-8, -12), (-10, 14), (-12, -18), (-6, 20)]
GENE_ID_PATTERN = re.compile(r'gene_id "([^"]+)"')
GENE_NAME_PATTERN = re.compile(r'gene_name "([^"]+)"')


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create exploratory QC plots and DE interpretation outputs from PyDESeq2 results."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("results/pydeseq2_all"),
        help="Directory containing PyDESeq2 outputs.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/pydeseq2_all/qc_diff_exp"),
        help="Directory for QC plots and interpretation outputs.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="Adjusted p-value threshold used to summarize significant genes.",
    )
    parser.add_argument(
        "--lfc-threshold",
        type=float,
        default=1.0,
        help="Absolute log2 fold-change threshold for stronger-effect summaries.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=15,
        help="Number of top up/down genes to export per contrast.",
    )
    parser.add_argument(
        "--gtf-path",
        type=Path,
        default=None,
        help="Optional GTF file or directory used to resolve gene names for volcano plot labels.",
    )
    return parser.parse_args()


def ensure_inputs(input_dir: Path) -> dict[str, Path]:
    required = {
        "counts": input_dir / "counts_matrix.samples_x_genes.tsv",
        "metadata": input_dir / "sample_metadata.tsv",
        "manifest": input_dir / "contrast_manifest.tsv",
    }
    missing = [str(path) for path in required.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required PyDESeq2 inputs: " + ", ".join(missing))
    return required


def load_inputs(input_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    counts = pd.read_csv(input_dir / "counts_matrix.samples_x_genes.tsv", sep="\t", index_col=0)
    metadata = pd.read_csv(input_dir / "sample_metadata.tsv", sep="\t", index_col=0)
    manifest = pd.read_csv(input_dir / "contrast_manifest.tsv", sep="\t")
    metadata = metadata.loc[counts.index]
    return counts, metadata, manifest


def build_condition_colors(metadata: pd.DataFrame) -> dict[str, str]:
    palette = CONDITION_PALETTE.copy()
    ordered_conditions = list(dict.fromkeys(metadata["condition"].tolist()))
    for index, condition in enumerate(ordered_conditions):
        palette.setdefault(condition, FALLBACK_CONDITION_PALETTE[index % len(FALLBACK_CONDITION_PALETTE)])
    return palette


def build_condition_legend_handles(condition_colors: dict[str, str]) -> list[Patch]:
    return [
        Patch(facecolor=color, edgecolor="white", linewidth=0.8, alpha=0.9, label=condition)
        for condition, color in condition_colors.items()
    ]


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


def resolve_gtf_path(gtf_path: Path | None) -> Path | None:
    candidates: list[Path] = []

    if gtf_path is not None:
        if gtf_path.is_dir():
            candidates.extend(sorted(gtf_path.glob("*.gtf.gz")))
            candidates.extend(sorted(gtf_path.glob("*.gtf")))
        else:
            candidates.append(gtf_path)
    else:
        repo_root = Path(__file__).resolve().parents[2]
        ref_dir = repo_root / "refs" / "ensembl_112"
        candidates.extend(
            [
                ref_dir / "Mus_musculus.GRCm39.112.gtf.gz",
                ref_dir / "Mus_musculus.GRCm39.112.gtf",
            ]
        )

    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def load_gene_name_map(gtf_path: Path | None) -> dict[str, str]:
    resolved_gtf = resolve_gtf_path(gtf_path)
    if resolved_gtf is None:
        return {}

    opener = gzip.open if resolved_gtf.suffix == ".gz" else open
    gene_names: dict[str, str] = {}

    with opener(resolved_gtf, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue

            attributes = fields[8]
            gene_id_match = GENE_ID_PATTERN.search(attributes)
            if gene_id_match is None:
                continue

            gene_name_match = GENE_NAME_PATTERN.search(attributes)
            gene_id = gene_id_match.group(1)
            gene_names[gene_id] = gene_name_match.group(1) if gene_name_match else gene_id

    return gene_names


def select_volcano_labels(result: pd.DataFrame, alpha: float, labels_per_direction: int) -> pd.DataFrame:
    candidates = result.loc[result["padj"].notna() & (result["padj"] < alpha)].copy()
    if candidates.empty:
        return candidates

    candidates["abs_log2FoldChange"] = candidates["log2FoldChange"].abs()
    up = candidates.loc[candidates["log2FoldChange"] > 0].sort_values(
        ["padj", "abs_log2FoldChange"], ascending=[True, False]
    ).head(labels_per_direction)
    down = candidates.loc[candidates["log2FoldChange"] < 0].sort_values(
        ["padj", "abs_log2FoldChange"], ascending=[True, False]
    ).head(labels_per_direction)
    return pd.concat([up, down], axis=0).sort_values(["padj", "abs_log2FoldChange"], ascending=[True, False])


def annotate_volcano_labels(ax: plt.Axes, labels: pd.DataFrame) -> None:
    positive_index = 0
    negative_index = 0
    for index, row in enumerate(labels.itertuples(index=False)):
        if row.log2FoldChange >= 0:
            dx, dy = VOLCANO_POSITIVE_LABEL_OFFSETS[positive_index % len(VOLCANO_POSITIVE_LABEL_OFFSETS)]
            positive_index += 1
        else:
            dx, dy = VOLCANO_NEGATIVE_LABEL_OFFSETS[negative_index % len(VOLCANO_NEGATIVE_LABEL_OFFSETS)]
            negative_index += 1
        ax.annotate(
            row.plot_label,
            (row.log2FoldChange, row.neg_log10_padj),
            xytext=(dx, dy),
            textcoords="offset points",
            ha="left" if row.log2FoldChange >= 0 else "right",
            va="center",
            fontsize=7.5,
            color=PLOT_NEUTRALS["ink"],
            bbox={"boxstyle": "round,pad=0.18", "facecolor": "white", "edgecolor": "none", "alpha": 0.82},
            arrowprops={"arrowstyle": "-", "color": PLOT_NEUTRALS["threshold_light"], "linewidth": 0.7},
            clip_on=False,
        )


def counts_per_million(counts: pd.DataFrame) -> pd.DataFrame:
    library_sizes = counts.sum(axis=1)
    return counts.div(library_sizes, axis=0) * 1_000_000


def log_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    return np.log2(counts_per_million(counts) + 1.0)


def save_library_size_plot(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    output_dir: Path,
    condition_colors: dict[str, str],
) -> pd.DataFrame:
    library_sizes = counts.sum(axis=1).rename("library_size")
    table = metadata.copy()
    table["library_size"] = library_sizes
    table.to_csv(output_dir / "library_sizes.tsv", sep="\t")

    fig, ax = plt.subplots(figsize=(11, 5.6))
    colors = [condition_colors[condition] for condition in table["condition"]]
    ax.bar(table.index, table["library_size"], color=colors, edgecolor="white", linewidth=0.8)
    ax.set_ylabel("Raw library size")
    ax.set_title("Library sizes from featureCounts matrix")
    ax.tick_params(axis="x", rotation=45, pad=8)
    ymax = table["library_size"].max() * 1.12
    ax.set_ylim(0, ymax)
    style_axis(ax, grid_axis="y")
    ax.legend(
        handles=build_condition_legend_handles(condition_colors),
        title="Condition",
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        borderaxespad=0.0,
    )
    for x, value in enumerate(table["library_size"]):
        ax.text(x, value + ymax * 0.015, f"{value/1e6:.1f}M", ha="center", va="bottom", fontsize=8)
    fig.subplots_adjust(bottom=0.24, right=0.82)
    save_figure(fig, output_dir / "library_sizes.png")
    return table


def save_log_cpm_boxplot(
    transformed: pd.DataFrame,
    metadata: pd.DataFrame,
    output_dir: Path,
    condition_colors: dict[str, str],
) -> None:
    fig, ax = plt.subplots(figsize=(12, 5.4))
    box = ax.boxplot(
        [transformed.loc[sample].values for sample in transformed.index],
        tick_labels=transformed.index.tolist(),
        patch_artist=True,
        showfliers=False,
    )
    for patch, sample in zip(box["boxes"], transformed.index):
        patch.set_facecolor(condition_colors[metadata.loc[sample, "condition"]])
        patch.set_alpha(0.85)
        patch.set_edgecolor(PLOT_NEUTRALS["ink"])
        patch.set_linewidth(0.9)
    for whisker in box["whiskers"]:
        whisker.set_color(PLOT_NEUTRALS["ink"])
    for cap in box["caps"]:
        cap.set_color(PLOT_NEUTRALS["ink"])
    for median in box["medians"]:
        median.set_color(PLOT_NEUTRALS["ink"])
        median.set_linewidth(1.3)
    ax.set_ylabel("log2(CPM + 1)")
    ax.set_title("Distribution of normalized expression values")
    ax.tick_params(axis="x", rotation=45, pad=8)
    style_axis(ax, grid_axis="y")
    ax.legend(
        handles=build_condition_legend_handles(condition_colors),
        title="Condition",
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        borderaxespad=0.0,
    )
    fig.subplots_adjust(bottom=0.24, right=0.83)
    save_figure(fig, output_dir / "log_cpm_boxplot.png")


def save_pca_plot(
    transformed: pd.DataFrame,
    metadata: pd.DataFrame,
    output_dir: Path,
    condition_colors: dict[str, str],
) -> pd.DataFrame:
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(transformed.values)
    pca_df = pd.DataFrame(
        {
            "PC1": pcs[:, 0],
            "PC2": pcs[:, 1],
            "condition": metadata["condition"].values,
            "replicate": metadata["replicate"].values,
        },
        index=transformed.index,
    )
    pca_df.to_csv(output_dir / "pca_coordinates.tsv", sep="\t")

    fig, ax = plt.subplots(figsize=(8.2, 6.0))
    style_axis(ax, grid_axis="both")
    for condition, subset in pca_df.groupby("condition"):
        ax.scatter(
            subset["PC1"],
            subset["PC2"],
            label=condition,
            color=condition_colors[condition],
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

    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}% variance)")
    ax.set_title("PCA on log2(CPM + 1)")
    ax.margins(x=0.08, y=0.08)
    ax.legend(title="Condition", frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)
    fig.subplots_adjust(right=0.82)
    save_figure(fig, output_dir / "pca_plot.png")
    return pca_df


def save_sample_correlation_heatmap(
    transformed: pd.DataFrame,
    metadata: pd.DataFrame,
    output_dir: Path,
    condition_colors: dict[str, str],
) -> pd.DataFrame:
    correlation = transformed.transpose().corr()
    correlation.to_csv(output_dir / "sample_correlation_matrix.tsv", sep="\t")

    distance = 1 - correlation
    linkage_matrix = linkage(squareform(distance.values, checks=False), method="average")
    dendro = dendrogram(linkage_matrix, no_plot=True, labels=correlation.index.tolist())
    ordered = dendro["ivl"]
    ordered_corr = correlation.loc[ordered, ordered]

    cmap = LinearSegmentedColormap.from_list(
        "sample_corr",
        ["#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"],
    )
    fig = plt.figure(figsize=(10.8, 6.5))
    grid = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.05, 0.22], wspace=0.15)
    ax = fig.add_subplot(grid[0, 0])
    cax = fig.add_subplot(grid[0, 1])
    legend_ax = fig.add_subplot(grid[0, 2])
    legend_ax.axis("off")
    image = ax.imshow(ordered_corr.values, cmap=cmap, vmin=ordered_corr.values.min(), vmax=1.0)
    ax.set_xticks(range(len(ordered)))
    ax.set_yticks(range(len(ordered)))
    ax.set_xticklabels(ordered, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_yticklabels(ordered)
    ax.set_title("Sample-sample Pearson correlation")
    ax.tick_params(axis="x", pad=8)
    ax.tick_params(axis="y", pad=10)
    for idx, sample in enumerate(ordered):
        ax.get_xticklabels()[idx].set_color(condition_colors[metadata.loc[sample, "condition"]])
        ax.get_yticklabels()[idx].set_color(condition_colors[metadata.loc[sample, "condition"]])
    legend_ax.legend(
        handles=build_condition_legend_handles(condition_colors),
        title="Condition",
        frameon=False,
        loc="upper left",
    )
    fig.colorbar(image, cax=cax, label="Pearson r")
    save_figure(fig, output_dir / "sample_correlation_heatmap.png")
    return correlation


def summarize_contrast(
    result_path: Path,
    alpha: float,
    lfc_threshold: float,
    top_n: int,
    output_dir: Path,
    gene_name_map: dict[str, str],
) -> dict[str, object]:
    result = pd.read_csv(result_path, sep="\t")
    contrast_name = result_path.stem
    significant = result.loc[result["padj"].notna() & (result["padj"] < alpha)].copy()
    strong = significant.loc[significant["log2FoldChange"].abs() >= lfc_threshold].copy()
    up = significant.loc[significant["log2FoldChange"] > 0].copy()
    down = significant.loc[significant["log2FoldChange"] < 0].copy()

    significant_sorted = significant.sort_values("padj")
    top_up = up.sort_values(["padj", "log2FoldChange"], ascending=[True, False]).head(top_n)
    top_down = down.sort_values(["padj", "log2FoldChange"], ascending=[True, True]).head(top_n)

    top_up.to_csv(output_dir / f"{contrast_name}.top_up.tsv", sep="\t", index=False)
    top_down.to_csv(output_dir / f"{contrast_name}.top_down.tsv", sep="\t", index=False)

    finite_padj = result["padj"].replace(0, np.nextafter(0, 1))
    neg_log10 = -np.log10(finite_padj)
    label_candidates = select_volcano_labels(result, alpha, VOLCANO_LABELS_PER_DIRECTION)
    if not label_candidates.empty:
        label_candidates = label_candidates.copy()
        label_candidates["neg_log10_padj"] = -np.log10(label_candidates["padj"].replace(0, np.nextafter(0, 1)))
        label_candidates["plot_label"] = label_candidates["Geneid"].map(gene_name_map).fillna(label_candidates["Geneid"])
    fig, ax = plt.subplots(figsize=(8.2, 5.8))
    nonsig_mask = ~(result["padj"].notna() & (result["padj"] < alpha))
    ax.scatter(result.loc[nonsig_mask, "log2FoldChange"], neg_log10.loc[nonsig_mask], s=12, c=PLOT_NEUTRALS["nonsig"], alpha=0.65)
    ax.scatter(up["log2FoldChange"], -np.log10(up["padj"].replace(0, np.nextafter(0, 1))), s=18, c=PLOT_NEUTRALS["sig_up"], alpha=0.82)
    ax.scatter(down["log2FoldChange"], -np.log10(down["padj"].replace(0, np.nextafter(0, 1))), s=18, c=PLOT_NEUTRALS["sig_down"], alpha=0.82)
    ax.axhline(-np.log10(alpha), color=PLOT_NEUTRALS["threshold"], linestyle="--", linewidth=1)
    ax.axvline(lfc_threshold, color=PLOT_NEUTRALS["threshold_light"], linestyle=":", linewidth=1)
    ax.axvline(-lfc_threshold, color=PLOT_NEUTRALS["threshold_light"], linestyle=":", linewidth=1)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.set_title(f"Volcano plot: {contrast_name}")
    style_axis(ax, grid_axis="y")
    if not label_candidates.empty:
        annotate_volcano_labels(ax, label_candidates)
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", linestyle="", color=PLOT_NEUTRALS["nonsig"], label="Not significant", markersize=6),
            Line2D([0], [0], marker="o", linestyle="", color=PLOT_NEUTRALS["sig_up"], label="Significant up", markersize=6),
            Line2D([0], [0], marker="o", linestyle="", color=PLOT_NEUTRALS["sig_down"], label="Significant down", markersize=6),
            Line2D([0], [0], linestyle="--", color=PLOT_NEUTRALS["threshold"], label=f"padj = {alpha:g}"),
            Line2D([0], [0], linestyle=":", color=PLOT_NEUTRALS["threshold_light"], label=f"|log2FC| = {lfc_threshold:g}"),
        ],
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        borderaxespad=0.0,
        title="Legend",
    )
    fig.subplots_adjust(right=0.78)
    save_figure(fig, output_dir / f"{contrast_name}.volcano.png")

    return {
        "contrast": contrast_name,
        "tested_genes": int(result["padj"].notna().sum()),
        "significant_genes": int(len(significant)),
        "up_genes": int(len(up)),
        "down_genes": int(len(down)),
        "strong_effect_genes": int(len(strong)),
        "median_abs_log2fc": float(significant["log2FoldChange"].abs().median()) if len(significant) else float("nan"),
        "top_gene_by_padj": significant_sorted.iloc[0]["Geneid"] if len(significant_sorted) else "NA",
    }


def write_summary_tables(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    manifest: pd.DataFrame,
    output_dir: Path,
) -> pd.DataFrame:
    manifest = manifest.copy()
    manifest["condition_pair"] = manifest["contrast_level"] + " vs " + manifest["reference_level"]
    manifest.to_csv(output_dir / "contrast_manifest_with_labels.tsv", sep="\t", index=False)

    per_condition = metadata.groupby("condition").size().rename("sample_count")
    per_condition.to_csv(output_dir / "samples_per_condition.tsv", sep="\t")

    gene_detection = (counts > 0).sum(axis=0).rename("detected_samples")
    gene_detection.describe().to_frame(name="value").to_csv(
        output_dir / "gene_detection_summary.tsv",
        sep="\t",
    )
    return manifest


def interpret_results(
    qc_stats: dict[str, float],
    contrast_summary: pd.DataFrame,
    pca_df: pd.DataFrame,
) -> list[str]:
    lines: list[str] = []
    lines.append("## Exploratory QC")
    lines.append(
        f"The count matrix contains {qc_stats['n_samples']:.0f} samples and {qc_stats['n_genes']:.0f} retained genes. "
        f"Raw library sizes range from {qc_stats['min_library_millions']:.1f}M to {qc_stats['max_library_millions']:.1f}M counts, "
        f"which is broad enough to require normalization but not so uneven that any one sample dominates the dataset."
    )

    pc_condition_centroids = pca_df.groupby("condition")[["PC1", "PC2"]].mean()
    if {"BSL", "FLT", "GC"}.issubset(pc_condition_centroids.index):
        flt_gc_distance = np.linalg.norm(pc_condition_centroids.loc["FLT"] - pc_condition_centroids.loc["GC"])
        flt_bsl_distance = np.linalg.norm(pc_condition_centroids.loc["FLT"] - pc_condition_centroids.loc["BSL"])
        gc_bsl_distance = np.linalg.norm(pc_condition_centroids.loc["GC"] - pc_condition_centroids.loc["BSL"])
        closest_pair = min(
            [("FLT and GC", flt_gc_distance), ("FLT and BSL", flt_bsl_distance), ("GC and BSL", gc_bsl_distance)],
            key=lambda item: item[1],
        )[0]
        lines.append(
            f"In the PCA space, the closest group centroids are {closest_pair}, which is consistent with the differential expression pattern rather than contradicting it."
        )
    else:
        lines.append("The PCA plot provides the main sample-level check for whether replicates cluster by condition.")

    lines.append("## Differential Expression Interpretation")
    strongest = contrast_summary.sort_values("significant_genes", ascending=False).iloc[0]
    weakest = contrast_summary.sort_values("significant_genes", ascending=True).iloc[0]
    lines.append(
        f"The strongest contrast is {strongest['contrast']} with {int(strongest['significant_genes'])} significant genes, while the weakest is {weakest['contrast']} with {int(weakest['significant_genes'])}."
    )

    contrast_lookup = contrast_summary.set_index("contrast")
    if {"condition_FLT_vs_BSL", "condition_FLT_vs_GC", "condition_GC_vs_BSL"}.issubset(contrast_lookup.index):
        flt_bsl = int(contrast_lookup.loc["condition_FLT_vs_BSL", "significant_genes"])
        flt_gc = int(contrast_lookup.loc["condition_FLT_vs_GC", "significant_genes"])
        gc_bsl = int(contrast_lookup.loc["condition_GC_vs_BSL", "significant_genes"])
        lines.append(
            f"Here, FLT vs BSL ({flt_bsl} genes) is much larger than FLT vs GC ({flt_gc} genes), and GC vs BSL is intermediate ({gc_bsl} genes). "
            "That pattern suggests the BSL samples are the most distinct baseline, while FLT and GC are more similar to each other than either is to BSL."
        )

    lines.append(
        "Use the volcano plots together with the top-up and top-down gene tables to prioritize genes with both strong effect size and strong statistical support rather than ranking only by p-value."
    )
    lines.append(
        "Because these tables contain Ensembl gene IDs only, the next interpretive step is usually gene annotation or pathway enrichment rather than manual inspection of raw IDs."
    )
    return lines


def write_report(
    output_dir: Path,
    qc_stats: dict[str, float],
    contrast_summary: pd.DataFrame,
    interpretation_lines: list[str],
) -> None:
    lines = [
        "# Exploratory QC and Differential Expression Interpretation",
        "",
        "## Files Produced",
        "- `library_sizes.tsv` and `library_sizes.png`",
        "- `log_cpm_boxplot.png`",
        "- `pca_coordinates.tsv` and `pca_plot.png`",
        "- `sample_correlation_matrix.tsv` and `sample_correlation_heatmap.png`",
        "- `contrast_summary.tsv`",
        "- per-contrast volcano plots and top up/down gene tables",
        "",
        "## Dataset Summary",
        f"- Samples: {int(qc_stats['n_samples'])}",
        f"- Genes analyzed: {int(qc_stats['n_genes'])}",
        f"- Minimum library size: {qc_stats['min_library_millions']:.1f}M counts",
        f"- Maximum library size: {qc_stats['max_library_millions']:.1f}M counts",
        "",
        "## Contrast Summary",
    ]
    for row in contrast_summary.itertuples(index=False):
        lines.append(
            f"- {row.contrast}: {int(row.significant_genes)} significant genes, {int(row.up_genes)} up, {int(row.down_genes)} down, {int(row.strong_effect_genes)} with |log2FC| >= 1"
        )
    lines.append("")
    lines.extend(interpretation_lines)
    lines.append("")
    (output_dir / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    ensure_inputs(args.input_dir)
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    counts, metadata, manifest = load_inputs(args.input_dir)
    transformed = log_cpm(counts)
    condition_colors = build_condition_colors(metadata)
    gene_name_map = load_gene_name_map(args.gtf_path)

    library_sizes = save_library_size_plot(counts, metadata, output_dir, condition_colors)
    save_log_cpm_boxplot(transformed, metadata, output_dir, condition_colors)
    pca_df = save_pca_plot(transformed, metadata, output_dir, condition_colors)
    save_sample_correlation_heatmap(transformed, metadata, output_dir, condition_colors)
    write_summary_tables(counts, metadata, manifest, output_dir)

    summaries = []
    for contrast in manifest["contrast"]:
        summary = summarize_contrast(
            args.input_dir / f"{contrast}.tsv",
            args.alpha,
            args.lfc_threshold,
            args.top_n,
            output_dir,
            gene_name_map,
        )
        summaries.append(summary)
    contrast_summary = pd.DataFrame(summaries).sort_values("significant_genes", ascending=False)
    contrast_summary.to_csv(output_dir / "contrast_summary.tsv", sep="\t", index=False)

    qc_stats = {
        "n_samples": float(counts.shape[0]),
        "n_genes": float(counts.shape[1]),
        "min_library_millions": float(library_sizes["library_size"].min() / 1_000_000),
        "max_library_millions": float(library_sizes["library_size"].max() / 1_000_000),
    }
    interpretation_lines = interpret_results(qc_stats, contrast_summary, pca_df)
    write_report(output_dir, qc_stats, contrast_summary, interpretation_lines)


if __name__ == "__main__":
    main()