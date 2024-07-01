#!/bin/bash
datum=$(date +"%Y%m%d")
touch $datum"_ISYmu1_Twist4_cutadaptlog.txt"
for  filename in ../fastq_twist4/ISYmu1/*R1_001.fastq.gz; do
	shortname="${filename:0:-16}"
	#ISYmu1:
	cutadapt -j 0 -g TCATGAGAGTTTCAA --discard-untrimmed -o $shortname"_5trim.fastq.gz" $filename >> $datum"_ISYmu1_Twist4_cutadaptlog.txt"
	cutadapt -j 0 -a TTTTTTCTTCCT -m 20 -M 20 --discard-untrimmed -o $shortname"_ISYmu1_Spacer.fastq.gz" $shortname"_5trim.fastq.gz" >> $datum"_ISYmu1_Twist4_cutadaptlog.txt"

done
for  filename in ../fastq_twist4/ISYmu1/*R2_001.fastq.gz; do
	shortname="${filename:0:-16}"
	cutadapt -j 0 -g AAGTACAAGTGGTAGA --discard-untrimmed -o $shortname"_3trim.fastq.gz" $filename >> $datum"_ISYmu1_Twist4_cutadaptlog.txt"
	cutadapt -j 0 -a CCATGCCGAAGC --discard-untrimmed -o $shortname"_ISYmu1_target.fastq.gz" $shortname"_3trim.fastq.gz" >> $datum"_ISYmu1_Twist4_cutadaptlog.txt"
done
