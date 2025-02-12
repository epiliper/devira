#!/bin/bash
# rasusa reads --num 1000000 -o PiV2-089908_S8_L001_R1.fastq.gz -o PiV2-089908_S8_L001_R2.fastq.gz PiV2-089908_S8_L001_R1_clean.fastq.gz PiV2-089908_S8_L001_R2_clean.fastq.gz -s 100

# Run rasusa to downsample reads
# Get read count of original files
r1_count=$(gzcat PiV2-089908_S8_L001_R1.fastq.gz | wc -l)
r1_count=$(($r1_count / 4))

r2_count=$(gzcat PiV2-089908_S8_L001_R2.fastq.gz | wc -l)
r2_count=$(($r2_count / 4))

# Compute total read count and 1% of it
count=$(($r1_count + $r2_count))
count=$(($count / 2))

# Compute 1% of the total count using bc for floating point calculation
sub_number=$(echo "$count * 0.01" | bc | awk '{print int($1)}')

# Run rasusa again with the computed number
rasusa reads --num $sub_number -o sub_PiV2-089908_S8_L001_R1.fastq.gz -o sub_PiV2-089908_S8_L001_R2.fastq.gz PiV2-089908_S8_L001_R1.fastq.gz PiV2-089908_S8_L001_R2.fastq.gz -s 100
