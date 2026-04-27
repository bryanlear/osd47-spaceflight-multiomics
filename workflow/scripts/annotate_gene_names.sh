#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: workflow/scripts/annotate_gene_names.sh <qc_diff_exp_dir> <gtf_or_ref_dir>

Annotates every *.top_up.tsv and *.top_down.tsv file in the given QC/DE output
directory with a gene_name column derived from the Ensembl GTF

e.g.,:
  workflow/scripts/annotate_gene_names.sh results/pydeseq2_all/qc_diff_exp refs/ensembl_112
  workflow/scripts/annotate_gene_names.sh results/pydeseq2_all/qc_diff_exp refs/ensembl_112/Mus_musculus.GRCm39.112.gtf.gz
EOF
}

resolve_gtf() {
    local gtf_input="$1"

    if [[ -d "$gtf_input" ]]; then
        local matches=()
        while IFS= read -r match; do
            matches+=("$match")
        done < <(find "$gtf_input" -maxdepth 1 \( -name '*.gtf.gz' -o -name '*.gtf' \) | sort)

        if [[ ${#matches[@]} -eq 0 ]]; then
            echo "No .gtf or .gtf.gz file found in: $gtf_input" >&2
            exit 1
        fi

        printf '%s\n' "${matches[0]}"
        return
    fi

    if [[ ! -f "$gtf_input" ]]; then
        echo "GTF path does not exist: $gtf_input" >&2
        exit 1
    fi

    printf '%s\n' "$gtf_input"
}

build_map() {
    local gtf_path="$1"
    local map_path="$2"
    local reader=(cat)

    if [[ "$gtf_path" == *.gz ]]; then
        reader=(zcat)
    fi

    "${reader[@]}" "$gtf_path" | awk -F '\t' '
        $3 == "gene" {
            gene_id = ""
            gene_name = ""

            if (match($9, /gene_id "[^"]+"/)) {
                gene_id = substr($9, RSTART + 9, RLENGTH - 10)
            }
            if (match($9, /gene_name "[^"]+"/)) {
                gene_name = substr($9, RSTART + 11, RLENGTH - 12)
            }

            if (gene_id != "") {
                if (gene_name == "") {
                    gene_name = gene_id
                }
                print gene_id "\t" gene_name
            }
        }
    ' | sort -u > "$map_path"
}

annotate_table() {
    local input_tsv="$1"
    local map_path="$2"
    local output_tsv="${input_tsv%.tsv}.with_gene_name.tsv"

    awk -F '\t' -v OFS='\t' '
        NR == FNR {
            gene_name[$1] = $2
            next
        }
        FNR == 1 {
            printf "%s\tgene_name", $1
            for (i = 2; i <= NF; i++) {
                printf "\t%s", $i
            }
            printf "\n"
            next
        }
        {
            resolved_name = ($1 in gene_name ? gene_name[$1] : "NA")
            printf "%s\t%s", $1, resolved_name
            for (i = 2; i <= NF; i++) {
                printf "\t%s", $i
            }
            printf "\n"
        }
    ' "$map_path" "$input_tsv" > "$output_tsv"

    printf 'Wrote %s\n' "$output_tsv"
}

main() {
    if [[ ${1:-} == "-h" || ${1:-} == "--help" ]]; then
        usage
        exit 0
    fi

    if [[ $# -ne 2 ]]; then
        usage >&2
        exit 1
    fi

    local input_dir="$1"
    local gtf_path
    local map_path
    local found_any=0

    if [[ ! -d "$input_dir" ]]; then
        echo "Input directory does not exist: $input_dir" >&2
        exit 1
    fi

    gtf_path="$(resolve_gtf "$2")"
    map_path="$(mktemp)"
    trap "rm -f \"$map_path\"" EXIT

    build_map "$gtf_path" "$map_path"

    while IFS= read -r table_path; do
        found_any=1
        annotate_table "$table_path" "$map_path"
    done < <(find "$input_dir" -maxdepth 1 \( -name '*.top_up.tsv' -o -name '*.top_down.tsv' \) | sort)

    if [[ $found_any -eq 0 ]]; then
        echo "No *.top_up.tsv or *.top_down.tsv files found in: $input_dir" >&2
        exit 1
    fi
}

main "$@"