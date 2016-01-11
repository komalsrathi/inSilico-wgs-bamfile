# align inSilico fastq files to gr37.fasta to create inSilico bam file
bwa_mem.pl -o SAMPLENAME -r gr37.fasta -s SAMPLENAME -t 2 SAMPLENAME/*.fastq.gz
