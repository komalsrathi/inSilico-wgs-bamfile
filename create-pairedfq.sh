#!/bin/bash

# this script will create test_bam_1.fastq and test_bam_2.fastq read pair files

inSilicoReadpairs.pl \
    -r hg19.fasta \
    -rl 150 \
    -d 1 \
    -o test_bam
