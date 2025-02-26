for i in {1..4}; do bwa index ref_${i}.fna; done; 
for i in {1..4}; do bwa mem ref_${i}.fna ../fastq/clean_fastq/sample-0${i}_1.fastq ../fastq/clean_fastq/sample-0${i}_2.fastq > aligned_reads_sample${i}.sam; done
