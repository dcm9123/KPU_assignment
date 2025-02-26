#!/bin/bash
conda activate unicycler;
for i in {1..4}; do unicycler -1 aligned_reads_sample${i}_R1.fastq -2 aligned_reads_sample${i}_R2.fastq --mode normal --threads 8 -o unicycler_output_sample${i}; done;
