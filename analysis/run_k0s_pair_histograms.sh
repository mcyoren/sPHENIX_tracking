#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage:"
  echo "  $0 input_directory input_file_pattern output_directory output_root_filename [tree_name]"
  echo
  echo "Example:"
  echo "  $0 /path/to/files 'pairTree_*.root' output k0s_histograms.root pairTree"
}

if [[ $# -lt 4 || $# -gt 5 ]]; then
  usage
  exit 1
fi

INPUT_DIRECTORY="$1"
INPUT_FILE_PATTERN="$2"
OUTPUT_DIRECTORY="$3"
OUTPUT_ROOT_FILENAME="$4"
TREE_NAME="${5:-pairTree}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_MACRO="${SCRIPT_DIR}/MakeK0sPairHistograms.C"

if [[ ! -f "${ROOT_MACRO}" ]]; then
  echo "ERROR: macro not found: ${ROOT_MACRO}"
  exit 2
fi

if [[ ! -d "${INPUT_DIRECTORY}" ]]; then
  echo "ERROR: input directory does not exist: ${INPUT_DIRECTORY}"
  exit 3
fi

mkdir -p "${OUTPUT_DIRECTORY}"

escape_root_string() {
  local value="$1"
  value="${value//\\/\\\\}"
  value="${value//\"/\\\"}"
  printf '%s' "${value}"
}

INPUT_DIRECTORY_ESCAPED="$(escape_root_string "${INPUT_DIRECTORY}")"
INPUT_FILE_PATTERN_ESCAPED="$(escape_root_string "${INPUT_FILE_PATTERN}")"
OUTPUT_DIRECTORY_ESCAPED="$(escape_root_string "${OUTPUT_DIRECTORY}")"
OUTPUT_ROOT_FILENAME_ESCAPED="$(escape_root_string "${OUTPUT_ROOT_FILENAME}")"
TREE_NAME_ESCAPED="$(escape_root_string "${TREE_NAME}")"
ROOT_MACRO_ESCAPED="$(escape_root_string "${ROOT_MACRO}")"

echo "Input directory : ${INPUT_DIRECTORY}"
echo "Input pattern   : ${INPUT_FILE_PATTERN}"
echo "Output directory: ${OUTPUT_DIRECTORY}"
echo "Output file     : ${OUTPUT_ROOT_FILENAME}"
echo "Tree name       : ${TREE_NAME}"

root -l -b -q \
  "${ROOT_MACRO_ESCAPED}(\"${INPUT_DIRECTORY_ESCAPED}\",\"${INPUT_FILE_PATTERN_ESCAPED}\",\"${OUTPUT_DIRECTORY_ESCAPED}\",\"${OUTPUT_ROOT_FILENAME_ESCAPED}\",\"${TREE_NAME_ESCAPED}\")"

echo "Finished: ${OUTPUT_DIRECTORY}/${OUTPUT_ROOT_FILENAME}"
