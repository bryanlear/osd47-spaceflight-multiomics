#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/../.." && pwd)"

default_python="$(command -v python 2>/dev/null || true)"
if [[ -z "$default_python" ]]; then
    default_python="$(command -v python3 2>/dev/null || true)"
fi
if [[ -z "$default_python" && -x "$repo_root/.venv/bin/python" ]]; then
    default_python="$repo_root/.venv/bin/python"
fi

python_bin="${PYTHON_BIN:-$default_python}"
abundance_path="${ABUNDANCE_PATH:-$repo_root/proteomics/results/fragpipe_tmt10_ms3_run3/tmt-report/abundance_protein_MD.tsv}"
combined_path="${COMBINED_PATH:-$repo_root/proteomics/results/fragpipe_tmt10_ms3_run3/combined_protein.tsv}"
output_dir="${OUTPUT_DIR:-$script_dir/results/fragpipe_tmt10_ms3_run3_protein_de}"
min_number_psm="${MIN_NUMBER_PSM:-2}"
min_max_pep_prob="${MIN_MAX_PEP_PROB:-0.99}"
min_combined_spectral_count="${MIN_COMBINED_SPECTRAL_COUNT:-0}"
min_reps_per_condition="${MIN_REPS_PER_CONDITION:-2}"
alpha="${ALPHA:-0.05}"

if [[ -z "$python_bin" || ! -x "$python_bin" ]]; then
    echo "Python interpreter not found at $python_bin" >&2
    exit 1
fi

echo "Using Python: $python_bin"

if ! "$python_bin" -c "import matplotlib, numpy, pandas, scipy, statsmodels" >/dev/null 2>&1; then
    echo "Missing Python packages. Install them with:" >&2
    echo "  $python_bin -m pip install -r $script_dir/requirements.protein_de.txt" >&2
    exit 1
fi

"$python_bin" "$script_dir/prepare_protein_de_inputs.py" \
    --abundance "$abundance_path" \
    --combined "$combined_path" \
    --output-dir "$output_dir"

"$python_bin" "$script_dir/fit_protein_linear_models.py" \
    --annotations "$output_dir/protein_annotations.tsv" \
    --abundance-matrix "$output_dir/protein_abundance_matrix.tsv" \
    --sample-metadata "$output_dir/sample_metadata.tsv" \
    --output-dir "$output_dir" \
    --min-number-psm "$min_number_psm" \
    --min-max-pep-prob "$min_max_pep_prob" \
    --min-combined-spectral-count "$min_combined_spectral_count" \
    --min-reps-per-condition "$min_reps_per_condition" \
    --alpha "$alpha" \
    "$@"

"$python_bin" "$script_dir/plot_protein_linear_model_results.py" \
    --input-dir "$output_dir" \
    --output-dir "$output_dir" \
    --min-number-psm "$min_number_psm" \
    --min-max-pep-prob "$min_max_pep_prob" \
    --alpha "$alpha"