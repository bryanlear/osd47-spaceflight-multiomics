#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
default_fragpipe_dir="${script_dir}/../fragpipe-24.0"
default_mzml_dir="${script_dir}/../mzML"
default_output_dir="${script_dir}/../results/fragpipe_tmt10_ms3"
default_experiment_name="OSD47_TMT10"
default_workflow_name="TMT10-MS3.workflow"
default_reference_tag="pool"
default_decoy_tag="rev_"
default_docker_image="eclipse-temurin:21-jre-jammy"

MZML_FILES=()
ANNOTATION_FILES=()

usage() {
  cat <<EOF
Usage: $(basename "$0") --fasta /absolute/path/to/database.fasta [options]

Run FragPipe headless on Linux using the repo-local FragPipe 24.0 install for the
OSD-47 single-plex TMT10 SPS-MS3 dataset.

Required:
  --fasta PATH              FragPipe-ready FASTA database to search against.
                            The file must already contain decoys with the
                            ${default_decoy_tag} prefix.

Optional:
  --fragpipe-dir DIR        FragPipe 24.0 installation directory.
                            Default: ${default_fragpipe_dir}
  --tools-dir DIR           Folder FragPipe should scan recursively for MSFragger
                            and IonQuant jars.
                            Default: <fragpipe-dir>
  --mzml-dir DIR            Directory containing the mzML fractions.
                            Default: ${default_mzml_dir}
  --output-dir DIR          Output directory for FragPipe results.
                            Default: ${default_output_dir}
  --experiment-name NAME    Manifest experiment name for this plex.
                            Default: ${default_experiment_name}
  --workflow-name NAME      Built-in FragPipe workflow to copy from the local install.
                            Default: ${default_workflow_name}
  --reference-tag TAG       Sample-name tag for the pooled reference channel.
                            Default: ${default_reference_tag}
  --threads N               Pass --threads N to FragPipe.
  --ram GB                  Pass --ram GB to FragPipe.
  --docker                  Run the Linux FragPipe workflow inside Docker.
                            Useful on macOS when no Linux machine is available.
  --docker-image IMAGE      Container image to use with --docker.
                            Default: ${default_docker_image}
  --dry-run                 Generate files and validate the run without starting a search.
  -h, --help                Show this help message.

The mzML directory must contain exactly one *annotation.txt file and all
fractions from the same TMT plex in the same folder.

The selected tools directory must contain at least one MSFragger jar, one
IonQuant jar, and one diaTracer jar somewhere under it.

If you downloaded a raw Ensembl or UniProt protein FASTA, prepare it first with
Philosopher so decoys are added before running this wrapper.
EOF
}

absolute_path() {
  local target="$1"
  local parent

  if [[ -d "${target}" ]]; then
    (
      cd -- "${target}"
      pwd -P
    )
    return 0
  fi

  parent="$(dirname -- "${target}")"
  (
    cd -- "${parent}"
    printf '%s/%s\n' "$(pwd -P)" "$(basename -- "${target}")"
  )
}

collect_mzml_files() {
  local mzml_dir="$1"
  local mzml_path

  MZML_FILES=()

  while IFS= read -r mzml_path; do
    [[ -n "${mzml_path}" ]] || continue
    MZML_FILES+=("${mzml_path}")
  done < <(find "${mzml_dir}" -maxdepth 1 -type f -iname '*.mzml' | sort)

  if [[ ${#MZML_FILES[@]} -eq 0 ]]; then
    echo "No .mzML files found in: ${mzml_dir}" >&2
    return 1
  fi
}

collect_annotation_files() {
  local mzml_dir="$1"
  local annotation_path

  ANNOTATION_FILES=()

  while IFS= read -r annotation_path; do
    [[ -n "${annotation_path}" ]] || continue
    ANNOTATION_FILES+=("${annotation_path}")
  done < <(find "${mzml_dir}" -maxdepth 1 -type f -name '*annotation.txt' | sort)

  if [[ ${#ANNOTATION_FILES[@]} -ne 1 ]]; then
    echo "FragPipe requires exactly one *annotation.txt file in ${mzml_dir}. Found ${#ANNOTATION_FILES[@]}." >&2
    return 1
  fi
}

write_manifest() {
  local manifest_path="$1"
  local experiment_name="$2"
  local mzml_path

  : > "${manifest_path}"

  for mzml_path in "${MZML_FILES[@]}"; do
    printf '%s\t%s\t\tDDA\n' "${mzml_path}" "${experiment_name}" >> "${manifest_path}"
  done
}

copy_workflow() {
  local fragpipe_dir="$1"
  local workflow_name="$2"
  local destination="$3"
  local workflow_source="${fragpipe_dir}/workflows/${workflow_name}"

  if [[ ! -f "${workflow_source}" ]]; then
    echo "Workflow file does not exist: ${workflow_source}" >&2
    return 1
  fi

  cp "${workflow_source}" "${destination}"
}

patch_workflow() {
  local source_workflow="$1"
  local patched_workflow="$2"
  local fasta_path="$3"
  local reference_tag="$4"
  local header_line=""
  local line

  {
    IFS= read -r header_line || true
    printf '%s\n' "${header_line}"
    printf 'database.db-path=%s\n' "${fasta_path}"

    while IFS= read -r line || [[ -n "${line}" ]]; do
      case "${line}" in
        database.db-path=*)
          continue
          ;;
        tmtintegrator.ref_tag=*)
          printf 'tmtintegrator.ref_tag=%s\n' "${reference_tag}"
          ;;
        *)
          printf '%s\n' "${line}"
          ;;
      esac
    done
  } < "${source_workflow}" > "${patched_workflow}"
}

require_search_tools() {
  local tools_dir="$1"
  local msfragger_match
  local ionquant_match
  local diatracer_match

  msfragger_match="$(find "${tools_dir}" -maxdepth 3 -type f -iname 'MSFragger*.jar' | head -n 1)"
  ionquant_match="$(find "${tools_dir}" -maxdepth 3 -type f -iname 'IonQuant*.jar' | head -n 1)"
  diatracer_match="$(find "${tools_dir}" -maxdepth 3 -type f \( -iname 'diaTracer*.jar' -o -iname 'diatracer*.jar' \) | head -n 1)"

  if [[ -z "${msfragger_match}" ]]; then
    echo "Missing MSFragger jar under ${tools_dir}. Download MSFragger into this tools directory before running FragPipe on Linux." >&2
    return 1
  fi

  if [[ -z "${ionquant_match}" ]]; then
    echo "Missing IonQuant jar under ${tools_dir}. Download IonQuant into this tools directory before running FragPipe on Linux." >&2
    return 1
  fi

  if [[ -z "${diatracer_match}" ]]; then
    echo "Missing diaTracer jar under ${tools_dir}. FragPipe 24.0 headless exits early without a valid diaTracer jar, even for this TMT10-MS3 workflow." >&2
    echo "Download an official diaTracer jar into this tools directory before running FragPipe." >&2
    return 1
  fi
}

require_fragpipe_ready_fasta() {
  local fasta_path="$1"
  local decoy_tag="$2"
  local fragpipe_dir="$3"
  local philosopher_bin="${fragpipe_dir}/tools/Philosopher/philosopher-v5.1.3-RC9"

  if ! grep -q "^>${decoy_tag}" "${fasta_path}"; then
    echo "FASTA does not look FragPipe-ready: no headers start with >${decoy_tag}." >&2
    echo "Prepare the raw FASTA with Philosopher before running this wrapper." >&2
    echo "Example on Linux:" >&2
    echo "  mkdir -p /absolute/path/to/database_workspace" >&2
    echo "  cd /absolute/path/to/database_workspace" >&2
    echo "  \"${philosopher_bin}\" workspace --init" >&2
    echo "  \"${philosopher_bin}\" database --custom \"${fasta_path}\" --contam" >&2
    return 1
  fi
}

collect_docker_mounts() {
  local repo_root="$1"
  shift

  local mount_parent
  local path
  local -a mounts=()
  local -a seen=()

  mounts+=("${repo_root}")
  seen+=("${repo_root}")

  for path in "$@"; do
    [[ -n "${path}" ]] || continue

    if [[ -d "${path}" ]]; then
      mount_parent="${path}"
    else
      mount_parent="$(dirname -- "${path}")"
    fi

    if [[ ! -d "${mount_parent}" ]]; then
      continue
    fi

    case "${path}" in
      "${repo_root}"|"${repo_root}"/*)
        continue
        ;;
    esac

    if [[ " ${seen[*]} " == *" ${mount_parent} "* ]]; then
      continue
    fi

    seen+=("${mount_parent}")
    mounts+=("${mount_parent}")
  done

  printf '%s\n' "${mounts[@]}"
}

run_inside_docker() {
  local repo_root="$1"
  local docker_image="$2"
  shift 2

  local section="mounts"
  local arg
  local mount_path
  local -a mount_sources=()
  local -a script_args=()
  local -a docker_cmd=(docker run --rm --platform linux/amd64)

  for arg in "$@"; do
    if [[ "${arg}" == "--" ]]; then
      section="args"
      continue
    fi

    if [[ "${section}" == "mounts" ]]; then
      mount_sources+=("${arg}")
    else
      script_args+=("${arg}")
    fi
  done

  while IFS= read -r mount_path; do
    [[ -n "${mount_path}" ]] || continue
    docker_cmd+=( -v "${mount_path}:${mount_path}" )
  done < <(collect_docker_mounts "${repo_root}" "${mount_sources[@]}")

  docker_cmd+=(
    -w "${repo_root}"
    "${docker_image}"
    bash
    "${script_dir}/$(basename -- "$0")"
    "${script_args[@]}"
  )

  echo "Running in Docker with image ${docker_image}" >&2
  exec "${docker_cmd[@]}"
}

main() {
  local fasta_path=""
  local fragpipe_dir="${default_fragpipe_dir}"
  local tools_dir=""
  local mzml_dir="${default_mzml_dir}"
  local output_dir="${default_output_dir}"
  local experiment_name="${default_experiment_name}"
  local workflow_name="${default_workflow_name}"
  local reference_tag="${default_reference_tag}"
  local threads=""
  local ram=""
  local use_docker=false
  local docker_image="${default_docker_image}"
  local dry_run=false
  local workflow_template
  local workflow_patched
  local manifest_path
  local workflow_stem
  local -a docker_script_args=()
  local fragpipe_args=()
  local fragpipe_bin
  local arg
  local repo_root

  while [[ $# -gt 0 ]]; do
    case "$1" in
      --fasta)
        fasta_path="$2"
        shift 2
        ;;
      --fragpipe-dir)
        fragpipe_dir="$2"
        shift 2
        ;;
      --tools-dir)
        tools_dir="$2"
        shift 2
        ;;
      --mzml-dir)
        mzml_dir="$2"
        shift 2
        ;;
      --output-dir)
        output_dir="$2"
        shift 2
        ;;
      --experiment-name)
        experiment_name="$2"
        shift 2
        ;;
      --workflow-name)
        workflow_name="$2"
        shift 2
        ;;
      --reference-tag)
        reference_tag="$2"
        shift 2
        ;;
      --threads)
        threads="$2"
        shift 2
        ;;
      --ram)
        ram="$2"
        shift 2
        ;;
      --docker)
        use_docker=true
        shift
        ;;
      --docker-image)
        docker_image="$2"
        shift 2
        ;;
      --dry-run)
        dry_run=true
        shift
        ;;
      -h|--help)
        usage
        return 0
        ;;
      *)
        echo "Unknown argument: $1" >&2
        usage >&2
        return 1
        ;;
    esac
  done

  repo_root="$(cd -- "${script_dir}/../.." && pwd -P)"

  if [[ -z "${fasta_path}" ]]; then
    echo "--fasta is required." >&2
    usage >&2
    return 1
  fi

  fasta_path="$(absolute_path "${fasta_path}")"

  if [[ -n "${fragpipe_dir}" ]]; then
    fragpipe_dir="$(absolute_path "${fragpipe_dir}")"
  fi

  if [[ -n "${tools_dir}" ]]; then
    tools_dir="$(absolute_path "${tools_dir}")"
  fi

  if [[ -n "${mzml_dir}" ]]; then
    mzml_dir="$(absolute_path "${mzml_dir}")"
  fi

  if [[ -n "${output_dir}" ]]; then
    output_dir="$(absolute_path "${output_dir}")"
  fi

  if [[ "${use_docker}" == true && "$(uname)" != "Linux" ]]; then
    if ! command -v docker >/dev/null 2>&1; then
      echo "--docker was requested but the docker command is not available." >&2
      return 1
    fi

    docker_script_args=(
      --fasta "${fasta_path}"
      --fragpipe-dir "${fragpipe_dir}"
      --mzml-dir "${mzml_dir}"
      --output-dir "${output_dir}"
      --experiment-name "${experiment_name}"
      --workflow-name "${workflow_name}"
      --reference-tag "${reference_tag}"
    )

    if [[ -n "${tools_dir}" ]]; then
      docker_script_args+=(--tools-dir "${tools_dir}")
    fi

    if [[ -n "${threads}" ]]; then
      docker_script_args+=(--threads "${threads}")
    fi

    if [[ -n "${ram}" ]]; then
      docker_script_args+=(--ram "${ram}")
    fi

    if [[ "${dry_run}" == true ]]; then
      docker_script_args+=(--dry-run)
    fi

    run_inside_docker "${repo_root}" "${docker_image}" \
      "${fasta_path}" "${fragpipe_dir}" "${mzml_dir}" "${output_dir}" "${tools_dir}" \
      -- \
      "${docker_script_args[@]}"
  fi

  if [[ ! -d "${fragpipe_dir}" ]]; then
    echo "FragPipe directory does not exist: ${fragpipe_dir}" >&2
    return 1
  fi

  if [[ ! -f "${fasta_path}" ]]; then
    echo "FASTA file does not exist: ${fasta_path}" >&2
    return 1
  fi

  if [[ ! -d "${mzml_dir}" ]]; then
    echo "mzML directory does not exist: ${mzml_dir}" >&2
    return 1
  fi

  mkdir -p "${output_dir}"

  if [[ -z "${tools_dir}" ]]; then
    tools_dir="${fragpipe_dir}"
  fi

  if [[ ! -d "${tools_dir}" ]]; then
    echo "Tools directory does not exist: ${tools_dir}" >&2
    return 1
  fi

  fragpipe_bin="${fragpipe_dir}/bin/fragpipe"

  require_fragpipe_ready_fasta "${fasta_path}" "${default_decoy_tag}" "${fragpipe_dir}"

  if [[ ! -x "${fragpipe_bin}" ]]; then
    echo "FragPipe launcher is not executable: ${fragpipe_bin}" >&2
    return 1
  fi

  if [[ "$(uname)" != "Linux" ]]; then
    echo "Warning: this wrapper targets Linux. Current platform is $(uname)." >&2
  fi

  workflow_stem="${workflow_name%.workflow}"
  workflow_template="${output_dir}/${workflow_stem}.template.workflow"
  workflow_patched="${output_dir}/${workflow_stem}.headless.workflow"
  manifest_path="${output_dir}/fragpipe_manifest.tsv"

  collect_mzml_files "${mzml_dir}"
  collect_annotation_files "${mzml_dir}"
  require_search_tools "${tools_dir}"
  write_manifest "${manifest_path}" "${experiment_name}"
  copy_workflow "${fragpipe_dir}" "${workflow_name}" "${workflow_template}"
  patch_workflow "${workflow_template}" "${workflow_patched}" "${fasta_path}" "${reference_tag}"

  fragpipe_args=(
    --headless
    --workflow
    "${workflow_patched}"
    --manifest
    "${manifest_path}"
    --workdir
    "${output_dir}"
    --config-tools-folder
    "${tools_dir}"
  )

  if [[ -n "${threads}" ]]; then
    fragpipe_args+=(--threads "${threads}")
  fi

  if [[ -n "${ram}" ]]; then
    fragpipe_args+=(--ram "${ram}")
  fi

  if [[ "${dry_run}" == true ]]; then
    fragpipe_args+=(--dry-run)
  fi

  echo "Using annotation file: ${ANNOTATION_FILES[0]}"
  echo "Wrote manifest: ${manifest_path}"
  echo "Wrote workflow: ${workflow_patched}"
  echo "Running FragPipe from ${fragpipe_bin}"

  for arg in "${fragpipe_args[@]}"; do
    printf ' %q' "${arg}"
  done
  printf '\n'

  "${fragpipe_bin}" "${fragpipe_args[@]}"
}

main "$@"