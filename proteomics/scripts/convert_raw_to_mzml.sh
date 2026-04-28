#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
default_input="${script_dir}/../extracted_raw"
default_output="${script_dir}/../mzML"
default_rosetta_env_name="thermo-raw-rosetta"

THERMO_RUNNER=()
SELECTED_CONVERTER=""

usage() {
  cat <<EOF
Usage: $(basename "$0") [input_path] [output_dir]

Convert Thermo .raw files to .mzML.

Defaults:
  input_path  ${default_input}
  output_dir  ${default_output}

Script prefers ThermoRawFileParser on macOS/Linux and falls back to
ProteoWizard msconvert if it is installed.

Optional environment variables:
  CONVERTER=thermorawfileparser|msconvert
  THERMO_RAW_ROSETTA_PREFIX=/path/to/thermo-raw-rosetta
  THERMO_RAW_FILE_PARSER_DLL=/path/to/ThermoRawFileParser.dll
EOF
}

resolve_rosetta_thermo_runner() {
  local rosetta_env_name="${THERMO_RAW_ROSETTA_ENV_NAME:-${default_rosetta_env_name}}"
  local candidate_prefixes=()
  local prefix

  if [[ -n "${THERMO_RAW_ROSETTA_PREFIX:-}" ]]; then
    candidate_prefixes+=("${THERMO_RAW_ROSETTA_PREFIX}")
  fi

  if [[ -n "${CONDA_PREFIX:-}" ]]; then
    candidate_prefixes+=("$(dirname -- "${CONDA_PREFIX}")/${rosetta_env_name}")
  fi

  if [[ -n "${CONDA_EXE:-}" ]]; then
    candidate_prefixes+=("${CONDA_EXE%/bin/conda}/envs/${rosetta_env_name}")
  fi

  candidate_prefixes+=(
    "${HOME}/anaconda3/envs/${rosetta_env_name}"
    "${HOME}/miniconda3/envs/${rosetta_env_name}"
    "${HOME}/miniforge3/envs/${rosetta_env_name}"
    "${HOME}/mambaforge/envs/${rosetta_env_name}"
  )

  for prefix in "${candidate_prefixes[@]}"; do
    [[ -n "${prefix}" ]] || continue

    if [[ -x "${prefix}/bin/ThermoRawFileParser" ]]; then
      THERMO_RUNNER=("${prefix}/bin/ThermoRawFileParser")
      return 0
    fi

    if [[ -x "${prefix}/bin/ThermoRawFileParser.sh" ]]; then
      THERMO_RUNNER=("${prefix}/bin/ThermoRawFileParser.sh")
      return 0
    fi
  done

  return 1
}

resolve_thermo_runner() {
  if resolve_rosetta_thermo_runner; then
    return 0
  fi

  if command -v ThermoRawFileParser >/dev/null 2>&1; then
    THERMO_RUNNER=(ThermoRawFileParser)
    return 0
  fi

  if command -v ThermoRawFileParser.sh >/dev/null 2>&1; then
    THERMO_RUNNER=(ThermoRawFileParser.sh)
    return 0
  fi

  if [[ -n "${THERMO_RAW_FILE_PARSER_DLL:-}" ]]; then
    if [[ ! -f "${THERMO_RAW_FILE_PARSER_DLL}" ]]; then
      echo "THERMO_RAW_FILE_PARSER_DLL does not exist: ${THERMO_RAW_FILE_PARSER_DLL}" >&2
      return 1
    fi

    if ! command -v dotnet >/dev/null 2>&1; then
      echo "dotnet is required to run THERMO_RAW_FILE_PARSER_DLL" >&2
      return 1
    fi

    THERMO_RUNNER=(dotnet "${THERMO_RAW_FILE_PARSER_DLL}")
    return 0
  fi

  return 1
}

resolve_converter() {
  case "${CONVERTER:-auto}" in
    auto)
      if resolve_thermo_runner; then
        SELECTED_CONVERTER="thermorawfileparser"
        return 0
      fi

      if command -v msconvert >/dev/null 2>&1; then
        SELECTED_CONVERTER="msconvert"
        return 0
      fi
      ;;
    thermorawfileparser)
      if resolve_thermo_runner; then
        SELECTED_CONVERTER="thermorawfileparser"
        return 0
      fi

      echo "CONVERTER=thermorawfileparser was requested, but ThermoRawFileParser was not found." >&2
      return 1
      ;;
    msconvert)
      if command -v msconvert >/dev/null 2>&1; then
        SELECTED_CONVERTER="msconvert"
        return 0
      fi

      echo "CONVERTER=msconvert was requested, but msconvert was not found." >&2
      return 1
      ;;
    *)
      echo "Unsupported CONVERTER value: ${CONVERTER}" >&2
      return 1
      ;;
  esac

  echo "No converter found. Install ThermoRawFileParser or msconvert, or set THERMO_RAW_FILE_PARSER_DLL." >&2
  return 1
}

collect_raw_files() {
  local input_path="$1"
  local candidate

  RAW_FILES=()

  if [[ ! -e "${input_path}" ]]; then
    echo "Input path does not exist: ${input_path}" >&2
    return 1
  fi

  if [[ "$(basename -- "${input_path}")" == *.[Rr][Aa][Ww] ]]; then
    RAW_FILES=("${input_path}")
    return 0
  fi

  if [[ ! -d "${input_path}" ]]; then
    echo "Input must be a .raw file or a directory containing .raw files: ${input_path}" >&2
    return 1
  fi

  while IFS= read -r -d '' candidate; do
    RAW_FILES+=("${candidate}")
  done < <(find "${input_path}" \( -type f -o -type d \) -iname '*.raw' -print0)

  if [[ ${#RAW_FILES[@]} -eq 0 ]]; then
    echo "No .raw files found under: ${input_path}" >&2
    return 1
  fi
}

convert_with_thermo_raw_file_parser() {
  local raw_path="$1"
  local output_dir="$2"

  "${THERMO_RUNNER[@]}" "-i=${raw_path}" "-o=${output_dir}" "-f=1" "-m=2" "-l=2"
}

convert_with_msconvert() {
  local raw_path="$1"
  local output_dir="$2"

  msconvert "${raw_path}" --mzML --zlib -o "${output_dir}"
}

main() {
  local input_path="${1:-${default_input}}"
  local output_dir="${2:-${default_output}}"
  local converter
  local raw_path
  local raw_name
  local raw_stem
  local output_file
  local converted_count=0
  local skipped_count=0

  case "${input_path}" in
    -h|--help)
      usage
      return 0
      ;;
  esac

  resolve_converter
  converter="${SELECTED_CONVERTER}"
  collect_raw_files "${input_path}"
  mkdir -p "${output_dir}"

  for raw_path in "${RAW_FILES[@]}"; do
    raw_name="$(basename -- "${raw_path}")"
    raw_stem="${raw_name%.[Rr][Aa][Ww]}"
    output_file="${output_dir}/${raw_stem}.mzML"

    if [[ -f "${output_file}" ]]; then
      echo "Skipping existing output: ${output_file}"
      skipped_count=$((skipped_count + 1))
      continue
    fi

    echo "Converting ${raw_name} with ${converter}"

    if [[ "${converter}" == "thermorawfileparser" ]]; then
      convert_with_thermo_raw_file_parser "${raw_path}" "${output_dir}"
    else
      convert_with_msconvert "${raw_path}" "${output_dir}"
    fi

    converted_count=$((converted_count + 1))
  done

  echo "Finished: converted ${converted_count} file(s), skipped ${skipped_count} existing file(s)."
}

RAW_FILES=()
main "$@"