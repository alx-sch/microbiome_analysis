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

run_mapping.sh (single FASTQ file):
```bash
#!/bin/bash

exp="SRX3198644"

mouse_index="1_index/mouse_index"
pathogen_index="1_index/pathogen_combined_index"
input_folder="2_exp_fastq"
output_folder="3_alignment/$exp"
stats_folder="4_stats/$exp"
input_reads="$input_folder/${exp}.fastq"

mkdir -p "$output_folder"
mkdir -p "$stats_folder"

# Align against the mouse genome
echo -e "-----\nStep 1: Aligning reads against the mouse genome...\n-----"
bwa mem -t 4 "$mouse_index" "$input_reads" | samtools view -b -f 4 -U "$output_folder/${exp}_mouse_mapped.bam" - > "$output_folder/${exp}_mouse_unmapped.bam"
echo -e "-----\nStep 1 completed.\n-----"

# Align unmapped-to-mouse reads against pathogens
echo -e "-----\nStep 2: Aligning unmapped-to-mouse reads against pathogens...\n-----"
samtools fastq -n -0 - "$output_folder/${exp}_mouse_unmapped.bam" | bwa mem -t 4 "$pathogen_index" - | samtools view -b -f 4 -U "$output_folder/${exp}_pathogen_mapped.bam" - > "$output_folder/${exp}_pathogen_unmapped.bam"
echo -e "-----\nStep 2 completed.\n-----"

# Sorting mapped pathogen reads and indexing
echo -e "-----\nStep 3: Sorting mapped pathogen reads and indexing...\n-----"
samtools sort -o "$output_folder/${exp}_pathogen_mapped_sorted.bam" "$output_folder/${exp}_pathogen_mapped.bam"
samtools index "$output_folder/${exp}_pathogen_mapped_sorted.bam"
echo -e "-----\nStep 3 completed.\n-----"

# Generating index statistics for mapped pathogen reads
echo -e "-----\nStep 4: Generating index statistics for mapped pathogen reads...\n-----"
samtools idxstats "$output_folder/${exp}_pathogen_mapped_sorted.bam" > "$stats_folder/${exp}_pathogen_mapped_idxstats.txt"
echo -e "-----\nStep 4 completed.\n-----\n"

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

### Display & save results

echo -e "Total reads: $reads_total (100%)\n-----" 
echo "Mouse: Mapped: $reads_mouse_mapped ($percentage_mouse_mapped%) | Unmapped: $reads_mouse_unmapped ($percentage_mouse_unmapped%)"
echo "Pathogen in non-mouse reads: Mapped: $reads_pathogen_mapped ($percentage_pathogen_mapped%) | Unmapped: $reads_pathogen_unmapped ($percentage_pathogen_unmapped%)"
echo "Total: Mapped: $reads_total_mapped ($percentage_total_mapped%) | Unmapped: $reads_pathogen_unmapped ($percentage_pathogen_unmapped%)"

echo -e "Total reads: $reads_total (100%)\n-----" > "$stats_folder/${exp}_mapping_results.txt"
echo "Mouse: Mapped: $reads_mouse_mapped ($percentage_mouse_mapped%) | Unmapped: $reads_mouse_unmapped ($percentage_mouse_unmapped%)" >> "$stats_folder/${exp}_mapping_results.txt"
echo "Pathogen in non-mouse reads: Mapped: $reads_pathogen_mapped ($percentage_pathogen_mapped%) | Unmapped: $reads_pathogen_unmapped ($percentage_pathogen_unmapped%)" >> "$stats_folder/${exp}_mapping_results.txt"
echo "Total: Mapped: $reads_total_mapped ($percentage_total_mapped%) | Unmapped: $reads_pathogen_unmapped ($percentage_pathogen_unmapped%)" >> "$stats_folder/${exp}_mapping_results.txt"
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

 ## Mapping to Pathogen Names

The idxstats report looks like this, containing multiple reference genomes (18) for the same pathogen (6):

```
NC_004917.1	1799146	0	0
NC_007795.1	2821361	0	0
NZ_KB944666.1	2806553	0	0
NZ_KB944667.1	4841	0	0
NZ_KB944668.1	58987	0	0
NZ_BBIX01000009.1	671026	0	0
NZ_BBIX01000008.1	854302	0	0
NZ_BBIX01000007.1	21295	0	0
NZ_BBIX01000006.1	814831	0	0
NZ_BBIX01000005.1	18018	0	0
NZ_BBIX01000004.1	14392	0	0
NZ_BBIX01000003.1	14977	0	0
NZ_BBIX01000002.1	11561	0	0
NZ_BBIX01000001.1	14144	0	0
NZ_CP033844.1	5864574	1	0
NZ_CP033845.1	6078	0	0
NZ_CP033846.1	8424	0	0
NZ_CP040863.1	2607636	0	0
*	0	0	0
```

pathogen_mapping.sh
```bash
#!/bin/bash

# Define the input file containing reference genomes
input_file="0_ref_genome/combined_pathogen_genome.fna"
# Define the output file
output_file="pathogen_mapping.txt"

# Declare a dictionary to store the mapping
declare -A reference_to_pathogen

# Extract reference genome identifiers and corresponding pathogen names
while IFS= read -r line; do
    # Extract reference genome identifier and pathogen name
    identifier=$(echo "$line" | awk -F' ' '{print $1}')
    pathogen=$(echo "$line" | awk -F' ' '{print $2 "_" $3}')
    
    # Store the mapping in the dictionary
    reference_to_pathogen["$identifier"]="$pathogen"
done < <(grep ">" "$input_file")

# Write the mapping to the output file
echo "Reference Genome Identifier	Pathogen Name" > "$output_file"
for key in "${!reference_to_pathogen[@]}"; do
    printf "%-30s %s\n" "$key" "${reference_to_pathogen[$key]}" >> "$output_file"
done
```
Dictionary: pathogen_mapping.txt
```
Reference Genome Identifier	Pathogen Name
>NZ_BBIX01000001.1             Rodentibacter_pneumotropicus
>NZ_BBIX01000007.1             Rodentibacter_pneumotropicus
>NZ_CP033846.1                 Klebsiella_oxytoca
>NZ_BBIX01000005.1             Rodentibacter_pneumotropicus
>NZ_CP040863.1                 Rodentibacter_heylii
>NZ_BBIX01000004.1             Rodentibacter_pneumotropicus
>NZ_KB944667.1                 Enterococcus_faecalis
>NZ_CP033844.1                 Klebsiella_oxytoca
>NZ_KB944666.1                 Enterococcus_faecalis
>NZ_BBIX01000008.1             Rodentibacter_pneumotropicus
>NC_007795.1                   Staphylococcus_aureus
>NZ_CP033845.1                 Klebsiella_oxytoca
>NZ_BBIX01000009.1             Rodentibacter_pneumotropicus
>NZ_KB944668.1                 Enterococcus_faecalis
>NZ_BBIX01000003.1             Rodentibacter_pneumotropicus
>NZ_BBIX01000006.1             Rodentibacter_pneumotropicus
>NZ_BBIX01000002.1             Rodentibacter_pneumotropicus
>NC_004917.1                   Helicobacter_hepaticus
```

 ## Results

```
Total Reads for subset_SRR26681942_1: 1116894 (100,00%)
------
Mouse: Mapped: 7790 (0,70%) | Unmapped: 1109104 (99,30%) 
Pathogen in non-mouse reads: Mapped: 23449 (2,10%) | Unmapped: 1085655 (97,20%)
Total: Mapped 31239 (2,80%) | Unmapped: 1085655 (97,20%)
```

- Depending on the experimental setup, skipping alignment against the mouse genome could significantly reduce processing time (it took ~80 min with around 1.1 M reads). The exclusion of this step doesn't result in a substantial increase in host-derived reads, accounting for only about 1% of the total reads in SRR26681942.
