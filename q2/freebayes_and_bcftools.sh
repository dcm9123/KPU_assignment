#!/bin/bash
for i in {1..4}; do freebayes -f reference_genes/gyrB_${i}.fasta fastq/clean_fastq/gyrB_aligned_sorted_${i}.bam > freebayes_output/gyrB_variants_${i}.vcf; done;
for i in {1..4}; do freebayes -f reference_genes/ppsA_${i}.fasta fastq/clean_fastq/ppsA_aligned_sorted_${i}.bam > freebayes_output/ppsA_variants_${i}.vcf; done;
for i in {1..4}; do bcftools filter -i 'QUAL>20 & INFO/DP>10' freebayes_output/gyrB_variants_${i}.vcf > freebayes_output/gyrB_variants_filtered_${i}.vcf; done;
for i in {1..4}; do bcftools filter -i 'QUAL>20 & INFO/DP>10' freebayes_output/ppsA_variants_${i}.vcf > freebayes_output/ppsA_variants_filtered_${i}.vcf; done;
