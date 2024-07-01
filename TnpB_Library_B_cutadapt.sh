#!/bin/bash
datum=$(date +"%Y%m%d")
touch $datum"_TnpB_cutadaptlog.txt"
for  filename in *R1_001.fastq.gz; do
	shortname="${filename:0:-21}"
	cutadapt -j 0 -g tcgtgtgaggttcaa --discard-untrimmed -o "Output/"$shortname"_5trim.fastq.gz" $filename >> $datum"_TnpB_cutadaptlog.txt"
	cutadapt -j 0 -a ggccggcatggt --discard-untrimmed -o "Output/"$shortname"_Spacer.fastq.gz" "Output/"$shortname"_5trim.fastq.gz" >> $datum"_TnpB_cutadaptlog.txt"
done
for  filename in *R2_001.fastq.gz; do
	shortname="${filename:0:-21}"
	cutadapt -j 0 -g GCTGAAGTACAAGTGGTAGACAGTC --discard-untrimmed -o "Output/"$shortname"_3trim.fastq.gz" $filename >> $datum"_TnpB_cutadaptlog.txt"
done
