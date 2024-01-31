nano # microbiome_analysis

- Download genomes:
  - Mouse: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/
  - Pathogens:
    - Helicobacter hepaticus: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000007905.1/
    - Staphylococcus aureus: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000013425.1/
    - Enterococcus faecalis: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000393015.1/
    - Rodentibacter pneumotropicus: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000730685.1/
    - Klebsiella oxytoca: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_003812925.1/
    - Rodentibacter heylii: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_010587025.1/
  
  - `datasets download genome accession GCF_000007905.1,GCF_000393015.1,GCF_003812925.1,GCF_000013425.1,GCF_000730685.1,GCF_010587025.1 --include genome`
-> .fna files

  - combine pathogen data: `cat GCF_000007905.1_ASM790v1_genomic.fna
GCF_000393015.1_EntefaecT5V1_genomic.fna GCF_003812925.1_ASM381292v1_genomic.fna
GCF_000013425.1_ASM1342v1_genomic.fna GCF_000730685.1ASM73068v1_genomic.fna
GCF_010587025.1_ASM1058702v1_genomic.fna> combined_pathogen_genomes.fna
`

- index:
   - `bwa index -p pathogen_combined_index combined_pathogen_genomes.fna`
   - `bwa index -p mouse_index GCF_000001635.27_GRCm39_genomic.fna` (~18 hours)

- download raw reads experimental data:
  - https://www.ebi.ac.uk/ena/browser/view/SRX3198644:  `fasterq-dump --outdir . --gzip --accession SRX3198644`
  - https://www.ncbi.nlm.nih.gov/sra/SRX22381836[accn]: `fastq-dump --outdir . --gzip --split-files SRR26681942`


pipeline:
1. align reads from exp (feces) to mouse genome
2. filter unmapped reads
3. align those reads to pathogens

```bash                                                               
#!/bin/bash

mouse_index="../1_index/mouse_index"
pathogen_index="../1_index/pathogen_combined_index"
input_reads="../2_exp_fastq/SRX3198644.fastq"

# Step 1: Align against the mouse genome
bwa mem -t 4 "$mouse_index" "$input_reads" > aligned_to_mouse.sam

# Step 2: Filter unmapped reads and convert to BAM
samtools view -b -f 4 aligned_to_mouse.sam > unmapped_to_mouse.bam

# Step 3: Align unmapped reads against the pathogens
bwa mem -t 4 "$pathogen_index" unmapped_to_mouse.bam | samtools view -b - > aligned_to_pathogen.bam

# Step 4: Generate index for the aligned_to_pathogen.bam file
samtools index aligned_to_pathogen.bam

# Step 5: Generate index statistics for the final BAM file
samtools idxstats aligned_to_pathogen.bam > aligned_to_pathogen_idxstats.txt
```


```
#!/bin/bash

mouse_index="../1_index/mouse_index"
pathogen_index="../1_index/pathogen_combined_index"
input_dir="../2_exp_fastq/"
output_dir="../3_alignment/"
stats_dir="../4_stats/"

# Check if there are any FASTQ files in the input directory
if [ -z "$(ls -A "$input_dir"/*.fastq)" ]; then
    echo "No FASTQ files found in $input_dir. Exiting."
    exit 1
fi

# Create directories if they don't exist
mkdir -p "$output_dir"
mkdir -p "$stats_dir"

for fastq_file in "$input_dir"/*.fastq; do
    # Extract the file name without extension
    file_name=$(basename -- "$fastq_file")
    file_name_no_ext="${file_name%.*}"

    # Create a subdirectory for each sample in the output and stats directories
    sample_output_dir="$output_dir$file_name_no_ext/"
    sample_stats_dir="$stats_dir$file_name_no_ext/"
    mkdir -p "$sample_output_dir"
    mkdir -p "$sample_stats_dir"

    # Step 1: Align against the mouse genome
    bwa mem -t 4 "$mouse_index" "$fastq_file" > "$sample_output_dir$file_name_no_ext"_aligned_to_mouse.sam

    # Step 2: Filter unmapped reads and convert to BAM
    samtools view -b -f 4 "$sample_output_dir$file_name_no_ext"_aligned_to_mouse.sam > "$sample_output_dir$file_name_no_ext"_unmapped_to_mouse.bam

    # Check if there are unmapped reads
    if [ -s "$sample_output_dir$file_name_no_ext"_unmapped_to_mouse.bam ]; then
        # Step 3: Align unmapped reads against the pathogens
        bwa mem -t 4 "$pathogen_index" "$sample_output_dir$file_name_no_ext"_unmapped_to_mouse.bam | samtools view -b - > "$sample_output_dir$file_name_no_ext"_aligned_to_pathogen.bam
        echo "Alignment against pathogens for $file_name completed in $sample_output_dir."

        # Step 4: Create index file for the aligned pathogen BAM
        samtools index "$sample_output_dir$file_name_no_ext"_aligned_to_pathogen.bam
        echo "Index file created for $file_name."

        # Step 5: Get index statistics for the aligned pathogen BAM file
        samtools idxstats "$sample_output_dir$file_name_no_ext"_aligned_to_pathogen.bam > "$sample_stats_dir$file_name_no_ext"_pathogen_idxstats.txt
        echo "Index statistics for $file_name written to $sample_stats_dir."
    else
        echo "No unmapped reads to pathogens for $file_name. Skipping pathogen alignment, index creation, and index stats steps."
    fi
done

echo "All samples processed."
```

#### BWA (Burrows-Wheeler Aligner):
Purpose:
- BWA is a software package used for aligning short DNA sequences against a large reference genome.
- It employs the Burrows-Wheeler Transform (BWT) algorithm to efficiently align short reads to a reference sequence.

Usage in the Script:
- In the script, BWA is used to align the experimental short DNA sequences (reads) against two different reference genomes: the mouse genome and a combined index for various pathogens.
- The `bwa mem` command is specifically used for DNA sequence alignment using the mem algorithm, suitable for short reads.
  
#### Samtools:
Purpose:
- Samtools is a suite of programs for interacting with high-throughput sequencing data in SAM (Sequence Alignment/Map) and BAM (Binary Alignment/Map) formats.
- It provides various utilities for manipulating, viewing, and analyzing sequence alignment data.
- 
Usage in the Script:
- In the script, Samtools is primarily used for two tasks:
  1. Filtering Unmapped Reads:
    - samtools view is used to filter reads based on specific criteria, such as selecting only unmapped reads ( `-f 4 ` flag).
    - The filtered reads are then redirected to a new BAM file (unmapped_to_mouse.bam).
      
#### Overall Workflow:
1. Alignment to Mouse Genome:
    - BWA is used to align short DNA reads from the experiment to the mouse genome, producing a SAM file (aligned_to_mouse.sam).

2. Filtering Unmapped Reads:
    - Samtools is employed to filter out unmapped reads from the SAM file, resulting in a new BAM file (unmapped_to_mouse.bam).

3. Conversion to BAM Format:
    - Samtools is used to explicitly convert the BAM file containing unmapped reads to BAM format.


4. Alignment to Pathogens:
   - BWA is used again to align the unmapped reads (now in BAM format) to a combined index for various pathogens, generating a BAM file (aligned_to_pathogen.bam).
Both BWA and Samtools play crucial roles in the processing and analysis of high-throughput sequencing data, especially in the context of DNA sequence alignment and manipulation.

- 

create index file:
`samtools index ./4_alignment/aligned_to_pathogen.bam`

then: `samtools idxstats ./4_alignment/aligned_to_pathogen.bam`


