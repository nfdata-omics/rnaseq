#!/bin/bash

# Default value for the working directory
DEFAULT_WORK_DIR="."

# Help message
if [ "$1" == "--help" ] || [ "$1" == "-h" ]; then
  echo "Usage: $0 [working_directory]"
  echo
  echo "Description:"
  echo "  This script generates an overall table with logFCs, raw pvalues and adjusted pvalues for all the tested comparisons, merging all the 'deseq2_toptable.*' files in the subfolders of the specified directory."
  echo "  If no directory is provided, the current directory is used by default."
  echo
  echo "Parameters:"
  echo "  working_directory (optional) - Path to the parent directory (usually dea_Set_Of_Samples) containing contrast subfolders (e.g. dea_Group-A-B) in which the 'deseq2_toptable.*' files to evaluate are stored."
  echo
  echo "Output:"
  echo "  An overall merged table named 'IPA_<working_directory>.txt' is created in the provided directory. This table is supposed to be used as input for the IPA analysis."
  exit 0
fi

# Check if a directory is provided as an argument, otherwise use the default
WORK_DIR="${1:-$DEFAULT_WORK_DIR}"

# Navigate to the provided working directory
cd "$WORK_DIR" || { echo "Error: Cannot navigate to directory $WORK_DIR"; exit 1; }


for i in dea_*/deseq2_toptable.*txt; do k=$(basename $i); j=${k%.txt}; z=${j#deseq2_toptable.}; cat <(cat $i | awk -v OFS='\t' -v z="$z" 'NR==1 {print z"."$1,z"."$3,z"."$4,z"."$5}') <(cat $i | awk -v OFS='\t' 'NR>1 {print $1,$3,$4,$5}'  | sort -k1,1) | cut -f 2- > ${z}_tmp; done;
cat <(cat $(ls dea_*/deseq2_toptable.*txt | head -n 1) | awk 'NR==1' | cut -f 1) <(cat $(ls dea_*/deseq2_toptable.*txt | head -n 1) | awk 'NR>1' | cut -f 1 | sort -k1,1) > genes;
paste genes *_tmp > IPA_$(basename $(pwd)).txt;
rm *_tmp genes;
