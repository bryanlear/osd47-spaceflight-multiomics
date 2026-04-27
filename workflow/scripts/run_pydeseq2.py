from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a PyDESeq2 contrast from a featureCounts gene count table."
    )
    parser.add_argument(
        "--counts",
        type=Path,
        default=Path("results/counts/osd47_featureCounts.txt"),
        help="Path to the featureCounts output table.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config/config.yaml"),
        help="Path to the pipeline config used to define sample metadata.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/pydeseq2"),
        help="Directory for PyDESeq2 outputs.",
    )
    parser.add_argument(
        "--design-factor",
        default="condition",
        help="Column in the metadata to use as the design factor.",
    )
    parser.add_argument(
        "--reference-level",
        default="BSL",
        help="Reference level for the design factor.",
    )
    parser.add_argument(
        "--contrast-level",
        default="FLT",
        help="Non-reference level to test against the reference level.",
    )
    parser.add_argument(
        "--all-pairwise",
        action="store_true",
        help=(
            "Run all pairwise contrasts across observed design-factor levels. "
            "Requested contrast level is treated as the primary group and is "
            "compared against every other group first."
        ),
    )
    parser.add_argument(
        "--min-count",
        type=int,
        default=10,
        help="Minimum count threshold for retaining a gene in a sample.",
    )
    parser.add_argument(
        "--min-samples",
        type=int,
        default=3,
        help="Minimum number of samples meeting --min-count for a gene to be kept.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="Adjusted p-value threshold for the significant-gene table.",
    )
    parser.add_argument(
        "--n-cpus",
        type=int,
        default=8,
        help="Number of CPU workers to give PyDESeq2.",
    )
    return parser.parse_args()


def load_config(config_path: Path) -> dict:
    with config_path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def sample_name_from_featurecounts_column(column_name: str) -> str:
    bam_name = Path(column_name).name
    if not bam_name.endswith(".sorted.bam"):
        raise ValueError(f"Unexpected featureCounts sample column: {column_name}")
    return bam_name.removesuffix(".sorted.bam")


def load_featurecounts_table(counts_path: Path) -> pd.DataFrame:
    raw_counts = pd.read_csv(counts_path, sep="\t", comment="#")
    count_columns = raw_counts.columns[6:]
    counts_by_gene = raw_counts.loc[:, ["Geneid", *count_columns]].copy()
    counts_by_gene = counts_by_gene.set_index("Geneid")
    counts_by_gene = counts_by_gene.rename(columns=sample_name_from_featurecounts_column)
    counts_by_gene = counts_by_gene.apply(pd.to_numeric)
    return counts_by_gene.transpose()


def build_metadata(config: dict, design_factor: str, sample_names: list[str]) -> pd.DataFrame:
    metadata = pd.DataFrame(
        [
            {
                "sample": archive["id"],
                "condition": archive["condition"],
                "replicate": archive["replicate"],
            }
            for archive in config["archives"]
        ]
    ).set_index("sample")

    missing_samples = sorted(set(sample_names) - set(metadata.index))
    if missing_samples:
        raise ValueError(
            "Samples were present in the featureCounts table but missing from config.yaml: "
            + ", ".join(missing_samples)
        )

    metadata = metadata.loc[sample_names].copy()
    if design_factor not in metadata.columns:
        raise ValueError(f"Design factor '{design_factor}' is not present in sample metadata.")
    return metadata


def order_factor_levels(
    metadata: pd.DataFrame,
    design_factor: str,
    reference_level: str,
    contrast_level: str,
) -> tuple[pd.DataFrame, list[str]]:
    observed_levels = list(pd.unique(metadata[design_factor]))
    if reference_level not in observed_levels:
        raise ValueError(
            f"Reference level '{reference_level}' is not present in {design_factor}: {observed_levels}"
        )
    if contrast_level not in observed_levels:
        raise ValueError(
            f"Contrast level '{contrast_level}' is not present in {design_factor}: {observed_levels}"
        )

    ordered_levels = [reference_level] + [
        level for level in observed_levels if level != reference_level
    ]
    metadata[design_factor] = pd.Categorical(
        metadata[design_factor],
        categories=ordered_levels,
        ordered=True,
    )
    return metadata, ordered_levels


def filter_low_count_genes(
    counts: pd.DataFrame,
    min_count: int,
    min_samples: int,
) -> pd.DataFrame:
    keep_mask = (counts >= min_count).sum(axis=0) >= min_samples
    filtered = counts.loc[:, keep_mask]
    if filtered.empty:
        raise ValueError(
            "All genes were filtered out. Relax --min-count or --min-samples for this dataset."
        )
    return filtered


def build_dds(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    design_factor: str,
    n_cpus: int,
) -> DeseqDataSet:
    design_formula = f"~{design_factor}"
    try:
        return DeseqDataSet(
            counts=counts,
            metadata=metadata,
            design=design_formula,
            refit_cooks=True,
            n_cpus=n_cpus,
        )
    except TypeError:
        return DeseqDataSet(
            counts=counts,
            metadata=metadata,
            design_factors=design_factor,
            refit_cooks=True,
            n_cpus=n_cpus,
        )


def build_contrasts(
    levels: list[str],
    reference_level: str,
    contrast_level: str,
    all_pairwise: bool,
) -> list[tuple[str, str]]:
    if not all_pairwise:
        return [(contrast_level, reference_level)]

    if len(levels) < 2:
        raise ValueError("At least two levels are required to run a contrast.")

    contrasts: list[tuple[str, str]] = []
    primary_pairs = [
        (contrast_level, level) for level in levels if level != contrast_level
    ]
    contrasts.extend(primary_pairs)

    remaining_levels = [level for level in levels if level != contrast_level]
    for lower_index, left_level in enumerate(remaining_levels[:-1]):
        for right_level in remaining_levels[lower_index + 1 :]:
            contrasts.append((right_level, left_level))

    return contrasts


def run_contrast(
    dds: DeseqDataSet,
    design_factor: str,
    contrast_level: str,
    reference_level: str,
    alpha: float,
    n_cpus: int,
    output_dir: Path,
) -> tuple[str, int]:
    stats = DeseqStats(
        dds,
        contrast=[design_factor, contrast_level, reference_level],
        alpha=alpha,
        n_cpus=n_cpus,
    )
    stats.summary()

    results = stats.results_df.sort_values("padj", na_position="last")
    contrast_stem = f"{design_factor}_{contrast_level}_vs_{reference_level}"
    results_path = output_dir / f"{contrast_stem}.tsv"
    results.to_csv(results_path, sep="\t")

    significant = results.loc[results["padj"].notna() & (results["padj"] < alpha)]
    significant.to_csv(
        output_dir / f"{contrast_stem}.significant.tsv",
        sep="\t",
    )
    return contrast_stem, len(significant)


def main() -> None:
    args = parse_args()
    config = load_config(args.config)

    counts = load_featurecounts_table(args.counts)
    metadata = build_metadata(config, args.design_factor, counts.index.tolist())
    metadata, observed_levels = order_factor_levels(
        metadata,
        args.design_factor,
        args.reference_level,
        args.contrast_level,
    )
    counts = filter_low_count_genes(counts, args.min_count, args.min_samples)
    metadata = metadata.loc[counts.index]

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    metadata.to_csv(output_dir / "sample_metadata.tsv", sep="\t")
    counts.to_csv(output_dir / "counts_matrix.samples_x_genes.tsv", sep="\t")

    dds = build_dds(counts, metadata, args.design_factor, args.n_cpus)
    dds.deseq2()

    contrasts = build_contrasts(
        observed_levels,
        args.reference_level,
        args.contrast_level,
        args.all_pairwise,
    )

    contrast_records = []
    for contrast_level, reference_level in contrasts:
        contrast_stem, significant_count = run_contrast(
            dds,
            args.design_factor,
            contrast_level,
            reference_level,
            args.alpha,
            args.n_cpus,
            output_dir,
        )
        contrast_records.append(
            {
                "contrast": contrast_stem,
                "contrast_level": contrast_level,
                "reference_level": reference_level,
                "significant_gene_count": significant_count,
            }
        )

    pd.DataFrame(contrast_records).to_csv(
        output_dir / "contrast_manifest.tsv",
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()