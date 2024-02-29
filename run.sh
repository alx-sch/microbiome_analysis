#!/bin/bash

exp="F11_S14_TRIMMED"

R1="${exp}_R1"
R2="${exp}_R2"

mouse_index="1_index/mouse_index"
pathogen_index="1_index/pathogen_combined_index"
input_folder="2_exp_fastq"
output_folder="3_alignment/$exp"
stats_folder="4_stats/$exp"
input_reads1="$input_folder/${R1}.fastq.gz"
input_reads2="$input_folder/${R2}.fastq.gz"

mkdir -p "$output_folder"
mkdir -p "$stats_folder"

# Check if mapping file exists. If not, create mapping file
if [ ! -f "4_stats/pathogen_mapping.txt" ]; then
    bash "utils/pathogen_mapping.sh"
fi

# Execute quack FASTQ QC
bash utils/quack_qc.sh ${exp}

# Align against the mouse genome
echo -e "-----\n${exp} / Step 1: Aligning reads against the mouse genome...\n-----"
bwa mem -t 4 "$mouse_index" "$input_reads1" "$input_reads2" | samtools view -b -f 4 -U "$output_folder/${exp}_mouse_mapped.bam" - > "$output_folder/${exp}_mouse_unmapped.bam"
echo -e "-----\nStep 1 completed.\n-----"

# Align unmapped-to-mouse reads against pathogens
echo -e "-----\n${exp} / Step 2: Aligning unmapped-to-mouse reads against pathogens...\n-----"
samtools fastq -n -0 - "$output_folder/${exp}_mouse_unmapped.bam" | bwa mem -t 4 "$pathogen_index" - | samtools view -b -f 4 -U "$output_folder/${exp}_pathogen_mapped.bam" - > "$output_folder/${exp}_pathogen_unmapped.bam"
echo -e "-----\nStep 2 completed.\n-----"

# Sorting mapped pathogen reads and indexing
echo -e "-----\n${exp} / Step 3: Sorting mapped pathogen reads and indexing...\n-----"
samtools sort -o "$output_folder/${exp}_pathogen_mapped_sorted.bam" "$output_folder/${exp}_pathogen_mapped.bam"
samtools index "$output_folder/${exp}_pathogen_mapped_sorted.bam"
echo -e "-----\nStep 3 completed.\n-----"

# Generating index statistics for mapped pathogen reads
echo -e "-----\n${exp} / Step 4: Generating index statistics for mapped pathogen reads...\n-----"
samtools idxstats "$output_folder/${exp}_pathogen_mapped_sorted.bam" > "$stats_folder/${exp}_idxstats.txt"

####
####

### Get number of reads
# Get number of reads (not) mapped to mouse (excl. supplementary (flag 2048) and secondary mapped reads (flag 256)).
reads_mouse_mapped=$(samtools view -c -F 2304 "$output_folder/${exp}_mouse_mapped.bam")
reads_mouse_unmapped=$(samtools view -c "$output_folder/${exp}_mouse_unmapped.bam")

# Calculate total reads
reads_total=$(( $reads_mouse_mapped + $reads_mouse_unmapped ))

# Get number of non-mouse reads (not) mapped to pathogens (excl. supplementary (flag 2048) and secondary mapped reads (flag 256)).
reads_pathogen_mapped=$(samtools view -c -F 2304 "$output_folder/${exp}_pathogen_mapped.bam")
reads_pathogen_unmapped=$(samtools view -c "$output_folder/${exp}_pathogen_unmapped.bam")

# Mouse and pathogens
reads_total_mapped=$((reads_mouse_mapped + reads_pathogen_mapped))

### Calculate percentages
percentage_mouse_mapped=$(awk "BEGIN {printf \"%.2f\", ($reads_mouse_mapped/$reads_total)*100}")
percentage_mouse_unmapped=$(awk "BEGIN {printf \"%.2f\", ($reads_mouse_unmapped/$reads_total)*100}")
percentage_pathogen_mapped=$(awk "BEGIN {printf \"%.2f\", ($reads_pathogen_mapped/$reads_total)*100}")
percentage_pathogen_unmapped=$(awk "BEGIN {printf \"%.2f\", ($reads_pathogen_unmapped/$reads_total)*100}")
percentage_total_mapped=$(awk "BEGIN {printf \"%.2f\", (($reads_total_mapped/$reads_total)*100)}")

####
####

### Translate reference genome identifiers to pathogen names and sum up reads
# Associative array to store the mapping
declare -A mapping

# Read mapping file and populate the dictionary
while IFS= read line; do
    identifier=$(echo "$line" | awk '{print $1}')
    name=$(echo "$line" | awk '{$1=""; print $0}' | xargs)
    mapping["$identifier"]="$name"
done < "4_stats/pathogen_mapping.txt"

# Associative array to store the mapped reads for each pathogen
declare -A pathogen_reads

while IFS=$'\t' read -r identifier length mapped_reads unmapped_reads; do
    # Get the translated pathogen name from the dictionary 
    translated="${mapping[$identifier]}"
    # If the translation exists, use it; otherwise, keep the original identifier
    if [ -n "$translated" ]; then
        (( pathogen_reads["$translated"] += mapped_reads ))
    else
        pathogen_reads["$identifier"]="$mapped_reads"   
    fi
done < "$stats_folder/${exp}_idxstats.txt"

echo -e "-----\nStep 4 completed.\n-----\n"

pathogen_sorted=("Enterococcus_faecalis" "Helicobacter_hepaticus" "Klebsiella_oxytoca" "Rodentibacter_heylii" "Rodentibacter_pneumotropicus" "Staphylococcus_aureus" "*")

# Write 'translation' and summed up reads into new file
echo -e "Total reads for ${exp}: $reads_total (100%)\n-----" > "$stats_folder/${exp}_summary.txt"
echo -e "Mouse:\t\t\tMapped: $reads_mouse_mapped ($percentage_mouse_mapped%)\t| Unmapped: $reads_mouse_unmapped ($percentage_mouse_unmapped%)" >> "$stats_folder/${exp}_summary.txt"
echo -e "Pathogen in non-mouse:\tMapped: $reads_pathogen_mapped ($percentage_pathogen_mapped%)\t| Unmapped: $reads_pathogen_unmapped ($percentage_pathogen_unmapped%)" >> "$stats_folder/${exp}_summary.txt"
echo -e "Total:\t\t\tMapped: $reads_total_mapped ($percentage_total_mapped%)\t| Unmapped: $reads_pathogen_unmapped ($percentage_pathogen_unmapped%)" >> "$stats_folder/${exp}_summary.txt"
echo -e "\n###########\n"  >> "$stats_folder/${exp}_summary.txt"
echo -e "Pathogen name\t\tMapped reads\t(% pathogen | % non-mouse) \n-----" >> "$stats_folder/${exp}_summary.txt"

for pathogen in "${pathogen_sorted[@]}"; do
    percentage_mapped=$(awk "BEGIN {printf \"%.2f\", (${pathogen_reads[$pathogen]}/$reads_pathogen_mapped)*100}")
    percentage_unmapped=$(awk "BEGIN {printf \"%.2f\", (${pathogen_reads[$pathogen]}/$reads_mouse_unmapped)*100}")
    echo -e "$pathogen\t\t${pathogen_reads[$pathogen]}\t(${percentage_mapped}%\t| ${percentage_unmapped}%)" >> "$stats_folder/${exp}_summary.txt"
 
done
