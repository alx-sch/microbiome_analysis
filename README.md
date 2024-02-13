# Detecting Pathogen DAN in Sequencing Data 

## System Specification
- **OS:** Ubuntu 22.04.01
- **Virtualization Platform:** Virtualbox 6.5.0
- **Architecture:** x86_64 GNU/Linux
- **Disk Space:**  100 GB
- **Memory (RAM):**  5.3 GB

## Reference Genomes
- Host (Mouse C57BL/6J):
  - https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/
- Pathogens:
  1. Helicobacter hepaticus: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000007905.1/
  2. Staphylococcus aureus: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000013425.1/
  3. Enterococcus faecalis: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000393015.1/
  4. Rodentibacter pneumotropicus: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000730685.1/
  5. Klebsiella oxytoca: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_003812925.1/
  6. Rodentibacter heylii: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_010587025.1/

Downloaded .fna files using 'ncbi-datasets' tool:
  - `datasets download genome accession GCF_000007905.1,GCF_000393015.1,GCF_003812925.1,GCF_000013425.1,GCF_000730685.1,GCF_010587025.1 --include genome`
  - Combine pathogen data: `cat GCF_000007905.1_ASM790v1_genomic.fna
GCF_000393015.1_EntefaecT5V1_genomic.fna GCF_003812925.1_ASM381292v1_genomic.fna
GCF_000013425.1_ASM1342v1_genomic.fna GCF_000730685.1ASM73068v1_genomic.fna
GCF_010587025.1_ASM1058702v1_genomic.fna> combined_pathogen_genomes.fna
`
- Excursion FASTA Files:
  - FASTA (.fna) files have the following structure:
    - Lines starting with '>' indicate the beginning of a sequence record (incl unqiue identifier).
    - The actual sequence data follows this indicator line (can be multiple lines)
    ```bash
    >NC_004917.1 Helicobacter hepaticus ATCC 51449, complete sequence
    CATTAAACCAAGTATAAAATCTATAAATTATCTTTA...
    ```
- Indexing:
  - BWA indexing is the process of creating an efficient data structure from a reference genome, allowing the BWA algorithm to rapidly align short DNA sequences during high-throughput sequencing analysis.
   - `bwa index -p pathogen_combined_index combined_pathogen_genomes.fna`
   - `bwa index -p mouse_index GCF_000001635.27_GRCm39_genomic.fna` (~18 hours!)

## Experimental Data

Download DNA sequencing data from experiments, onbtained from feces; ideally containing full genomes (microbial and host components):
  - https://www.ebi.ac.uk/ena/browser/view/SRX3198644:  `fasterq-dump --outdir . --gzip SRX3198644`
    - Only 74.9 kB --> 220 reads (`wc -l SRX3198644.fastq | awk '{print $1/4}'`)

  - https://www.ncbi.nlm.nih.gov/sra/SRX22381836[accn]: `fastq-dump --outdir . --gzip --split-files SRR26681942`
    - Both files 10.0 GB each (!)  --> two files as data was generated by pair-end sequencing
    - -> 22.4M reads each set --> heavy on computation --> create subsets for a general overview

  - Excursion: FASTQ Files:
    - Contains sequence data and quality information, as produced by sequencing machines
    - Each entry conists of four lines:
        - Sequence Indentifier Line: Starts with **'@'** followed by unique identifier for the read, may contain additional information like sequencing instrumet, sample details, read length, etc.
        - Sequence Line: Containts actual nucleotide sequence of DNA/RNA read.
        - Quality Score Identifier Line: Start with **'+'** and usually mirrors the sequence identifier line.
        - Quality Score Line: Contains ASCII-encoded quality scores, representing the confidence/accuracy of each base call.
      ```bash
      @SRX3198644.1 1 length=45
      TTGTTGAACTGGCTCTTTTTCGCAATCCCGCTGTAAGTACTGTCT
      +SRX3198644.1 1 length=45
      AAAAAEEEEEEEEEEEEEEEEEEEAEEAEEEEEEEEEEEEEEEEE
      ```

  -  Generating subsets of ~1M reads each using 'seqtk':
      - `seqtk sample -s 42 SRR26681942_1.fastq 0.05 > subset_SRR26681942_1.fastq`
      - `seqtk sample -s 42 SRR26681942_2.fastq 0.05 > subset_SRR26681942_2.fastq`
      - '42' is a seed value (making random sampling reproducale), new subsets contain 5% of total reads.

## Aligning Seq Data to Genomes

#### Suggested Directory Structure

- alignments/
	- run_alignment_all.sh
	- 0_ref_genome/  
		- mouse_genome.fna   
		- pathogen_1_genome.fna
  		- pathogen_2_genome.fna
  		- ...
        - pathogen_combined_genome.fna 	 	
 	 - 1_index/
   		- mouse_index.amb
     	- mouse_index.ann
       	- ...
       	- pathogen_combined_index.amb
       	- pathogen_combined_index.ann
       	- ...
     - 2_exp_fastq/
       	  -    exp1.fastq
       	  -    exp2.fastq
       	  -    ...
 


#### Pipeline
- Aligns experimental sequencing data reads to the mouse genome and pathogenic indices.
- Generates comprehensive statistical reports.

The script iterates over a set of FASTQ files, aligns them to the mouse genome, filters unmapped reads, converts them to FASTQ format, aligns them to a combined pathogen index, and performs various post-processing steps. Finally, it calculates and reports statistics on the alignment results, such as the percentage of reads mapped to the mouse genome, pathogens, and those not mapped to either. The results and statistics are saved in the specified output and stats folders.

run_alignement_all.sh:
```bash
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
```

run_alignement_single.sh:
```bash
#!/bin/bash

exp="SRX3198644"

input_folder="2_exp_fastq"
output_folder="3_alignment"
stats_folder="4_stats"

mkdir -p "$output_folder"
mkdir -p "$stats_folder"

# Align against the mouse genome
echo -e "-----\nStep 1: Aligning against the mouse genome...\n-----"
bwa mem -t 4 "1_index/mouse_index" "$input_reads" > "3_alignment/${exp}/${exp}_aligned_to_mouse.sam"
echo -e "-----\nStep 1 completed.\n-----"

# Filter unmapped reads
echo -e "-----\nStep 2: Filtering unmapped reads...\n-----"
samtools view -b -f 4 "3_alignment/${exp}/${exp}_aligned_to_mouse.sam" > "3_alignment/${exp}/${exp}_unmapped_to_mouse.bam"
echo "Step 2 completed.\n"

# Convert to FASTQ
echo -e "-----\nStep 3: Conversion to FASTQ...\n-----"
samtools fastq -n -0 "3_alignment/${exp}/${exp}_unmapped_to_mouse.fastq" "3_alignment/${exp}/${exp}_unmapped_to_mouse.bam"
echo -e "-----\nStep 3 completed.\n-----"

# Align unmapped reads against the pathogens
echo -e "-----\nStep 4: Align against pathogens...\n-----"
bwa mem -t 4 "1_index/pathogen_combined_index" "3_alignment/${exp}/${exp}_unmapped_to_mouse.fastq" > "3_alignment/${exp}/${exp}_aligned_to_pathogen.sam"
echo -e "-----\nStep 4 completed.\n-----"

# Convert SAM to BAM
echo -e "-----\nStep 5: Converting SAM to BAM...\n-----"
samtools view -b -o "$3_alignment/${exp}/${exp}_aligned_to_pathogen.bam" "3_alignment/${exp}/${exp}_aligned_to_pathogen.sam"
echo -e "-----\nStep 5 completed.\n-----"

# Sort BAM file
echo -e "-----\nStep 6: Sorting BAM file...\n-----"
samtools sort -o "3_alignment/${exp}/${exp}_aligned_to_pathogen_sorted.bam" "3_alignment/${exp}/${exp}_aligned_to_pathogen.bam"
echo -e "-----\nStep 6 completed.\n-----"

# Index the sorted BAM file
echo -e "-----\nStep 7: Indexing the sorted BAM file...\n-----"
samtools index "3_alignment/${exp}/${exp}_aligned_to_pathogen_sorted.bam"
echo -e "-----\nStep 7 completed.\n-----"

# Generate index statistics for the final BAM file
echo -e "-----\nStep 8: Generating index statistics...\n-----"
samtools idxstats "3_alignment/${exp}/${exp}_aligned_to_pathogen_sorted.bam" > "4_stats/${exp}/${exp}_aligned_to_pathogen_idxstats.txt"
echo -e "-----\nStep 8 completed.\n-----"

# Calculate total number of reads (excl. supplementary (flag 2048) and secondary reads (flag 256)).
total_reads=$(samtools view -c -F 2304 "3_alignment/${exp}/${exp}_aligned_to_mouse.sam")

# Calculate number of reads not mapped to mouse
unmapped_to_mouse_reads=$(samtools view -c -f 4 "3_alignment/${exp}/${exp}_aligned_to_mouse.sam")

# Calculate number of reads mapped to pathogens (excl. supplementary (flag 2048) and secondary reads (flag 256) along with unmapped (flag 4))
mapped_to_pathogen_reads=$(samtools view -c -F 2308 "3_alignment/${exp}/${exp}_aligned_to_pathogen_sorted.bam")

# Calculate number of reads not mapped to mouse or pathogens
unmapped_to_either_reads=$((unmapped_to_mouse_reads - mapped_to_pathogen_reads))

# Calculate percentages
percentage_not_mapped_to_mouse=$(awk "BEGIN {printf \"%.2f\", (${unmapped_to_mouse_reads}/${total_reads})*100}")
percentage_mapped_to_pathogen=$(awk "BEGIN {printf \"%.2f\", (${mapped_to_pathogen_reads}/${total_reads})*100}")
percentage_not_mapped_to_either=$(awk "BEGIN {printf \"%.2f\", (${unmapped_to_either_reads}/${total_reads})*100}")

# Display results
echo "Total Reads: $total_reads (100,00%)"
echo "------"
echo "Reads not mapped to mouse: $unmapped_to_mouse_reads (${percentage_not_mapped_to_mouse}%)"
echo "Reads mapped to pathogens: $mapped_to_pathogen_reads (${percentage_mapped_to_pathogen}%)"
echo "Reads not mapped to mouse or pathogen: $unmapped_to_either_reads (${percentage_not_mapped_to_either}%)"

# Redirect results to text file
echo -e "Total Reads: $total_reads (100,00%)\n------\nReads not mapped to mouse: $unmapped_to_mouse_reads (${percentage_not_mapped_to_mouse}%)\nReads mapped to pathogens: $mapped_to_pathogen_reads (${percentage_mapped_to_pathogen}%)\nReads not mapped to mouse or pathogen: $unmapped_to_either_reads (${percentage_not_mapped_to_either}%)" > "4_stats/${exp}/${exp}_alignment_results.txt"
```

run_alignement_single.sh // IMPR:
```bash
#!/bin/bash

exp="SRX3198644"

input_folder="2_exp_fastq"
output_folder="3_alignment"
stats_folder="4_stats"
input_reads="${input_folder}/${exp}.fastq"

mkdir -p "$output_folder"
mkdir -p "$stats_folder"

# Align against the mouse genome
echo -e "-----\nStep 1: Aligning against the mouse genome...\n-----"
bwa mem -t 4 "1_index/mouse_index" "$input_reads" > "3_alignment/${exp}/${exp}_aligned_to_mouse.sam"
echo -e "-----\nStep 1 completed.\n-----"

# Filter unmapped reads
echo -e "-----\nStep 2: Filtering unmapped reads...\n-----"
samtools view -b -f 4 -o "3_alignment/${exp}/${exp}_unmapped_to_mouse.bam" -U "3_alignment/${exp}/${exp}_mapped_to_mouse.bam" "3_alignment/${exp}/${exp}_aligned_to_mouse.sam"
echo "Step 2 completed.\n"

# Convert to FASTQ
echo -e "-----\nStep 3: Conversion to FASTQ...\n-----"
samtools fastq -n -0 "3_alignment/${exp}/${exp}_unmapped_to_mouse.fastq" "3_alignment/${exp}/${exp}_unmapped_to_mouse.bam"
echo -e "-----\nStep 3 completed.\n-----"

# Align unmapped reads against the pathogens
echo -e "-----\nStep 4: Align against pathogens...\n-----"
bwa mem -t 4 "1_index/pathogen_combined_index" "3_alignment/${exp}/${exp}_unmapped_to_mouse.fastq" > "3_alignment/${exp}/${exp}_aligned_to_pathogen.sam"
echo -e "-----\nStep 4 completed.\n-----"

# Convert SAM to BAM
echo -e "-----\nStep 5: Converting SAM to BAM...\n-----"
samtools view -b -o "3_alignment/${exp}/${exp}_aligned_to_pathogen.bam" -U "3_alignment/${exp}/${exp}_unaligned_to_pathogen.bam" "3_alignment/${exp}/${exp}_aligned_to_pathogen.sam"
echo -e "-----\nStep 5 completed.\n-----"

# Sort BAM file
echo -e "-----\nStep 6: Sorting BAM file...\n-----"
samtools sort -o "3_alignment/${exp}/${exp}_aligned_to_pathogen_sorted.bam" "3_alignment/${exp}/${exp}_aligned_to_pathogen.bam"
echo -e "-----\nStep 6 completed.\n-----"

# Index the sorted BAM file
echo -e "-----\nStep 7: Indexing the sorted BAM file...\n-----"
samtools index "3_alignment/${exp}/${exp}_aligned_to_pathogen_sorted.bam"
echo -e "-----\nStep 7 completed.\n-----"

# Generate index statistics for the final BAM file
echo -e "-----\nStep 8: Generating index statistics...\n-----"
samtools idxstats "3_alignment/${exp}/${exp}_aligned_to_pathogen_sorted.bam" > "4_stats/${exp}/${exp}_aligned_to_pathogen_idxstats.txt"
echo -e "-----\nStep 8 completed.\n-----"

# Calculate total number of reads (excl. supplementary (flag 2048) and secondary reads (flag 256)).
total_reads=$(samtools view -c -F 2304 "3_alignment/${exp}/${exp}_aligned_to_mouse.sam")

# Calculate number of reads not mapped to mouse
unmapped_to_mouse_reads=$(samtools view -c -f 4 "3_alignment/${exp}/${exp}_aligned_to_mouse.sam")

# Calculate number of reads mapped to pathogens (excl. supplementary (flag 2048) and secondary reads (flag 256) along with unmapped (flag 4))
mapped_to_pathogen_reads=$(samtools view -c -F 2308 "3_alignment/${exp}/${exp}_aligned_to_pathogen_sorted.bam")

# Calculate number of reads not mapped to mouse or pathogens
unmapped_to_either_reads=$((unmapped_to_mouse_reads - mapped_to_pathogen_reads))

# Calculate percentages
percentage_not_mapped_to_mouse=$(awk "BEGIN {printf \"%.2f\", (${unmapped_to_mouse_reads}/${total_reads})*100}")
percentage_mapped_to_pathogen=$(awk "BEGIN {printf \"%.2f\", (${mapped_to_pathogen_reads}/${total_reads})*100}")
percentage_not_mapped_to_either=$(awk "BEGIN {printf \"%.2f\", (${unmapped_to_either_reads}/${total_reads})*100}")

## Display results
echo "Total Reads: $total_reads (100.00%)"
echo "------"
echo "Reads mapped to mouse: Y: $mapped_to_mouse_reads ($(printf "%.2f" $percentage_mapped_to_mouse)%) / N: $unmapped_to_mouse_reads ($(printf "%.2f" $percentage_not_mapped_to_mouse)%)"
echo "Reads mapped to pathogens: Y: $mapped_to_pathogen_reads ($(printf "%.2f" $percentage_mapped_to_pathogen)%) / N: $unmapped_to_either_reads ($(printf "%.2f" $percentage_not_mapped_to_either)%)"
echo "Reads mapped to neither: Y: 0 (0.00%) / N: $unmapped_to_either_reads (100.00%)"

# Redirect results to text file
echo -e "Total Reads: $total_reads (100.00%)\n------\nReads mapped to mouse: Y: $mapped_to_mouse_reads ($(printf "%.2f" $percentage_mapped_to_mouse)%) / N: $unmapped_to_mouse_reads ($(printf "%.2f" $percentage_not_mapped_to_mouse)%)\nReads mapped to pathogens: Y: $mapped_to_pathogen_reads ($(printf "%.2f" $percentage_mapped_to_pathogen)%) / N: $unmapped_to_either_reads ($(printf "%.2f" $percentage_not_mapped_to_either)%)\nReads mapped to neither: Y: 0 (0.00%) / N: $unmapped_to_either_reads (100.00%)" > "4_stats/${exp}/${exp}_alignment_results.txt"
```

#### Tools

**BWA (Burrows-Wheeler Aligner):**
- Purpose:
	- BWA is a software package used for aligning short DNA sequences against a large reference genome.
	- It employs the Burrows-Wheeler Transform (BWT) algorithm to efficiently align short reads to a reference sequence.
- Usage in script:
	- BWA is used to align the experimental short DNA sequences (reads) against two different reference genomes: the mouse genome and a combined index for various pathogens.
	- The `bwa mem` command is specifically used for DNA sequence alignment using the mem algorithm, suitable for short reads.
 	- The alignment results are stored in SAM (Sequence Alignment/Map) format files.
  
**Samtools:**
- Purpose:
	- Samtools is a suite of programs for interacting with high-throughput sequencing data in SAM (Sequence Alignment/Map) and BAM (Binary Alignment/Map) formats.
	- It provides various utilities for manipulating, viewing, and analyzing sequence alignment data.
- Usage in script:
   	- Filtering unmapped reads from the alignment results.
   	- Converting SAM to BAM format.
   	- Indexing the sorted BAM file.
	- Generating index statistics for the final BAM file.
	- Calculating the number of reads mapped to different categories (mouse, pathogens, neither).

 ## Results
```
Total Reads for SRX3198644: 220 (100,00%)
------
Reads not mapped to mouse: 193 (87,73%)
Reads mapped to pathogens: 1 (0,45%)
Reads not mapped to mouse or pathogen: 192 (87,27%)
```
```
Total Reads for subset_SRR26681942_1: 1116893 (100,00%)
------
Reads not mapped to mouse: 1109104 (99,30%)
Reads mapped to pathogens: 23449 (2,10%)
Reads not mapped to mouse or pathogen: 1085655 (97,20%)
```
```
Total Reads for subset_SRR26681942_2: 1116893 (100,00%)
------
Reads not mapped to mouse: 1105597 (98,99%)
Reads mapped to pathogens: 23281 (2,08%)
Reads not mapped to mouse or pathogen: 1082316 (96,90%)
```
- Depending on the experimental setup, skipping alignment against the mouse genome could significantly reduce processing time (it took ~80 min with around 1.1 M reads). The exclusion of this step doesn't result in a substantial increase in host-derived reads, accounting for only about 1% of the total reads in SRR26681942.
