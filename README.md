# microbiome_analysis


download raw reads experimental data:
 `fasterq-dump --outdir /path/to/output_directory --gzip --accession SRX3198644 `
 https://www.ebi.ac.uk/ena/browser/view/SRX3198644

pipeline:
1. align reads from exp (feces) to mouse genome
2. filter unmapped reads
3. align those reads to pathogens

```bash                                                               
#!/bin/bash

mouse_index="../1_index/mouse_index"
pathogen_index="../1_index/pathogen_combined_index"
input_reads="../3_exp_seq/SRX3198644.fastq"

# Step 1: Align against the mouse genome
bwa mem -t 4 "$mouse_index" "$input_reads" > aligned_to_mouse.sam

# Step 2: Filter unmapped reads
samtools view -b -f 4 aligned_to_mouse.sam > unmapped_to_mouse.bam
# Convert unmapped reads to BAM format explicitly
samtools view -b unmapped_to_mouse.bam > unmapped_to_mouse.bam

# Step 4: Align unmapped reads against the pathogens
bwa mem -t 4 "$pathogen_index" unmapped_to_mouse.bam | samtools view -b - > aligned_to_pathogen.bam
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
  2. Converting SAM to BAM:
    - `samtools view ` is again used, this time to convert a SAM file (unmapped reads) to BAM format.
    - The result is stored in the same BAM file (unmapped_to_mouse.bam), ensuring it is in BAM format for further processing.
      
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
