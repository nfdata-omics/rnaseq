#!/bin/bash

# Default value for the working directory
DEFAULT_WORK_DIR="."

# Help message
if [ "$1" == "--help" ] || [ "$1" == "-h" ]; then
  echo "Usage: $0 [working_directory]"
  echo
  echo "Description:"
  echo "  This script generates an overall summary file from a set of 'deseq2_summary.*' files in the specified directory."
  echo "  If no directory is provided, the current directory is used by default."
  echo
  echo "Parameters:"
  echo "  working_directory (optional) - Path to the directory containing the 'deseq2_summary.*' files to evaluate."
  echo
  echo "Output:"
  echo "  A summary file named '<working_directory>_summary' is created in the provided directory."
  exit 0
fi

# Check if a directory is provided as an argument, otherwise use the default
WORK_DIR="${1:-$DEFAULT_WORK_DIR}"

# Navigate to the provided working directory
cd "$WORK_DIR" || { echo "Error: Cannot navigate to directory $WORK_DIR"; exit 1; }

# Producing the overall dge summary for the considered set of samples
cat <(echo -e 'Comparison\tUpreg\tUpreg_percent\tDownreg\tDownreg_percent') <(for i in deseq2_summary.*; do j=${i%.txt}; z=${j#deseq2_summary.}; paste <(echo $z) <(cat $i | grep up | sed 's/.*: //' | tr -d ' ' | tr ',' '\t') <(cat $i | grep down | sed 's/.*: //' | tr -d ' ' | tr ',' '\t'); done) > $(basename $(pwd))_summary
