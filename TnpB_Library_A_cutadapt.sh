#!/bin/bash
datum=$(date +"%Y%m%d")
touch $datum"_TnpB_cutadaptlog.txt"
for  filename in *R1_001.fastq.gz; do
	shortname="${filename:8:-20}"
	cutadapt -j 0 -g TCGTGTGAGGTTCAA -o "Output/"$shortname"_5trim.fastq.gz" $filename >> $datum"_TnpB_cutadaptlog.txt"
	cutadapt -j 0 -a GGCCGGCATGGTCCC --discard-untrimmed -o "Output/"$shortname"_Spacer.fastq.gz" "Output/"$shortname"_5trim.fastq.gz" >> $datum"_TnpB_cutadaptlog.txt"
done
for  filename in *R2_001.fastq.gz; do
	shortname="${filename:8:-20}"
	cutadapt -j 0 -g AGAGCCATTTGTCTCGCTGA --discard-untrimmed -o "Output/"$shortname"_3trim.fastq.gz" $filename >> $datum"_TnpB_cutadaptlog.txt"
done
