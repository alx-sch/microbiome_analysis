#!/bin/bash

input_folder="2_exp_fastq"
output_folder="3_alignment"
stats_folder="4_stats"

mkdir -p "$output_folder"
mkdir -p "$stats_folder"

# Iterate over all fastq files in the input folder
for fastq_file in "$input_folder"/*.fastq; do

	# Extract the experiment name from the file name
    	exp=$(basename "$fastq_file" .fastq)
	input_reads="${input_folder}/${exp}.fastq"

	mkdir -p "$output_folder/${exp}"
	mkdir -p "$stats_folder/${exp}"

	# Align against the mouse genome
   	echo -e "-----\nStep 1: Aligning against the mouse genome for ${exp}...\n-----"
    	bwa mem -t 4 "1_index/mouse_index" "$input_reads" > "$output_folder/$exp/${exp}_aligned_to_mouse.sam"
    	echo -e "-----\nStep 1 completed.\n-----"

	# Filter unmapped reads
    	echo -e "-----\nStep 2: Filtering unmapped reads for ${exp}...\n-----"
    	samtools view -b -f 4 "$output_folder/$exp/${exp}_aligned_to_mouse.sam" > "$output_folder/$exp/${exp}_unmapped_to_mouse.bam"
    	echo "Step 2 completed.\n"

	# Convert to FASTQ
    	echo -e "-----\nStep 3: Conversion to FASTQ for ${exp}...\n-----"
    	samtools fastq -n -0 "$output_folder/$exp/${exp}_unmapped_to_mouse.fastq" "$output_folder/$exp/${exp}_unmapped_to_mouse.bam"
    	echo -e "-----\nStep 3 completed.\n-----"

 	# Align unmapped reads against the pathogens
    	echo -e "-----\nStep 4: Align against pathogens for ${exp}...\n-----"
    	bwa mem -t 4 "1_index/pathogen_combined_index" "$output_folder/$exp/${exp}_unmapped_to_mouse.fastq" > "$output_folder/$exp/${exp}_aligned_to_pathogen.sam"
    	echo -e "-----\nStep 4 completed.\n-----"

	# Convert SAM to BAM
    	echo -e "-----\nStep 5: Converting SAM to BAM for ${exp}...\n-----"
    	samtools view -b -o "$output_folder/$exp/${exp}_aligned_to_pathogen.bam" "$output_folder/$exp/${exp}_aligned_to_pathogen.sam"
    	echo -e "-----\nStep 5 completed.\n-----"

	# Sort BAM file
    	echo -e "-----\nStep 6: Sorting BAM file for ${exp}...\n-----"
    	samtools sort -o "$output_folder/$exp/${exp}_aligned_to_pathogen_sorted.bam" "$output_folder/$exp/${exp}_aligned_to_pathogen.bam"
    	echo -e "-----\nStep 6 completed.\n-----"

	# Index the sorted BAM file
    	echo -e "-----\nStep 7: Indexing the sorted BAM file for ${exp}...\n-----"
    	samtools index "$output_folder/$exp/${exp}_aligned_to_pathogen_sorted.bam"
    	echo -e "-----\nStep 7 completed.\n-----"

	# Generate index statistics for the final BAM file
    	echo -e "-----\nStep 8: Generating index statistics for ${exp}...\n-----"
    	samtools idxstats "$output_folder/$exp/${exp}_aligned_to_pathogen_sorted.bam" > "$stats_folder/$exp/${exp}_aligned_to_pathogen_idxstats.txt"
    	echo -e "-----\nStep 8 completed.\n-----"

	# Calculate total number of reads (excl. supplementary (flag 2048) and secondary reads (flag 256)).
    	total_reads=$(samtools view -c -F 2304 "$output_folder/$exp/${exp}_aligned_to_mouse.sam")

	# Calculate number of reads not mapped to mouse
    	unmapped_to_mouse_reads=$(samtools view -c -f 4 "$output_folder/$exp/${exp}_aligned_to_mouse.sam")

	# Calculate number of reads mapped to pathogens (excl. supplementary (flag 2048) and secondary reads (flag 256) along with unmapped (flag 4))
    	mapped_to_pathogen_reads=$(samtools view -c -F 2308 "$output_folder/$exp/${exp}_aligned_to_pathogen_sorted.bam")

	# Calculate number of reads not mapped to mouse or pathogens
    	unmapped_to_either_reads=$((unmapped_to_mouse_reads - mapped_to_pathogen_reads))

	# Calculate percentages
    	percentage_not_mapped_to_mouse=$(awk "BEGIN {printf \"%.2f\", (${unmapped_to_mouse_reads}/${total_reads})*100}")
    	percentage_mapped_to_pathogen=$(awk "BEGIN {printf \"%.2f\", (${mapped_to_pathogen_reads}/${total_reads})*100}")
    	percentage_not_mapped_to_either=$(awk "BEGIN {printf \"%.2f\", (${unmapped_to_either_reads}/${total_reads})*100}")

	# Display and save current report
    	report="Total Reads for ${exp}: $total_reads (100,00%)\n------\nReads not mapped to mouse: $unmapped_to_mouse_reads (${percentage_not_mapped_to_mouse}%)\nReads mapped to pathogens: $mapped_to_pathogen_reads (${percentage_mapped_to_pathogen}%)\nReads not mapped to mouse or pathogen: $unmapped_to_either_reads (${percentage_not_mapped_to_either}%)"
	echo -e "$report"
	echo -e "$report" > "$stats_folder/$exp/${exp}_alignment_results.txt"

	done
