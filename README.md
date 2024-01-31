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

mouse_index="1_index/mouse_index"
pathogen_index="1_index/pathogen_combined_index"
input_reads="3_exp_seq/SRX3198644.fastq"

# Step 1: Align against the mouse genome
bwa mem -t 4 "$mouse_index" "$input_reads" > aligned_to_mouse.sam

# Step 2: Filter unmapped reads
samtools view -b -f 4 aligned_to_mouse.sam > unmapped_to_mouse.bam
# Convert unmapped reads to BAM format explicitly
samtools view -b unmapped_to_mouse.bam > unmapped_to_mouse.bam

# Step 4: Align unmapped reads against the pathogens
bwa mem -t 4 "$pathogen_index" unmapped_to_mouse.bam | samtools view -b - > aligned_to_pathogen.bam
```

