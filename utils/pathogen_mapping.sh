#!/bin/bash

input_file="0_ref_genome/combined_pathogen_genome.fna"
output_file="4_stats/pathogen_mapping.txt"

mkdir -p "4_stats"

# Declare a dictionary to store the mapping
declare -A reference_to_pathogen

# Extract reference genome identifiers and corresponding pathogen names
while IFS= read -r line; do
    identifier=$(echo "$line" | awk -F' ' '{print $1}')
    pathogen=$(echo "$line" | awk -F' ' '{print $2 "_" $3}')
    
    # Store the mapping in the dictionary
    reference_to_pathogen["$identifier"]="$pathogen"
done < <(grep ">" "$input_file")

# Write the mapping to output file
echo "Reference Genome Identifier	Pathogen Name" > "$output_file"
echo -e "Mapping file created..."
for key in "${!reference_to_pathogen[@]}"; do
    cleaned_key="${key#>}"
    printf "%-30s %s\n" "$cleaned_key" "${reference_to_pathogen[$key]}" >> "$output_file"
done
