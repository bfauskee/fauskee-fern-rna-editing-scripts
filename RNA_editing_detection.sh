#!/bin/bash



#Updated RNA editiing site detection protocol from Edera and Sanchez-Puerta 2020

#First assess trimmed reads with FastQC. If per base sequence content is noisy at the front of the read use the following second trim step.
#Run Trimmomatic in SE mode on F and R reads seperately.

TAXON=""

module load Trimmomatic

java -jar $TRIMMOMATIC SE /work/bdf24/plastomes/RNA_editing/ophioglossales/reads/${TAXON}_RNA_PE_1.fastq /work/bdf24/plastomes/RNA_editing/ophioglossales/reads/${TAXON}_1_paired_GCtrim.fastq HEADCROP:13 CROP:72
java -jar $TRIMMOMATIC SE /work/bdf24/plastomes/RNA_editing/ophioglossales/reads/${TAXON}_RNA_PE_2.fastq /work/bdf24/plastomes/RNA_editing/ophioglossales/reads/${TAXON}_2_paired_GCtrim.fastq HEADCROP:13 CROP:72

#cd /work/bdf24/plastomes/RNA_editing
#parse discription lines to make human friendly
sed -i 's/.∗gene=\([^]]∗\)] .∗/>\1/g' ${TAXON}_CDS.fasta

#Build bowtie index
module load Bowtie2
bowtie2-build ${TAXON}_CDS.fasta ${TAXON}_CDS

#Run Bowtie2 to map RNA reads to CDS multi-fasta
bowtie2 -x ${TAXON}_CDS -1 /work/bdf24/plastomes/RNA_editing/ophioglossales/reads/${TAXON}_1_paired_GCtrim.fastq -2 /work/bdf24/plastomes/RNA_editing/ophioglossales/reads/${TAXON}_2_paired_GCtrim.fastq \
-S ${TAXON}.sam --local --no-unal --no-mixed --fr --nofw -p 64

#Now compress SAM file into BAM file
module load samtools
samtools view -S -b ${TAXON}.sam > ${TAXON}.bam

#Now sort the bam, and index the sorted bam.
samtools sort ${TAXON}.bam -o ${TAXON}.srt.bam
samtools index ${TAXON}.srt.bam

#Starting RNA editing site detection
#BAM readcount will extract read nucleotides alinged to each reference position
#RUN THIS IN BACKGROUND
module load bam-readcount
bam-readcount -f ${TAXON}_CDS.fasta ${TAXON}.srt.bam > ${TAXON}.rc 

#BAM-readcount does not include positions with no aligned reads so we will have to create this file from the first .rc file
#we create a blank BAM-readcount file from fasta file, this file has blanks for every single ref. position.

cat ${TAXON}_CDS.fasta | sed 's/^>\([^\n]\+\)/>\1!/g' | tr -d '\n' | tr '>' '!' | sed 's/^!//g' | awk 'BEGIN{K=""; P=0; T=":0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00"} {split($0, a, "!"); for (i=1; i<=length(a); i+=2) {split(a[i+1], seq,""); for (j=1; j<=length(seq); j++) printf "%s\t%d\t%s\t%d\t=%s\tA%s\tC%s\tG%s\tT%s\tN%s\n", a[i], j, seq[j], 0, T, T, T, T, T, T}}' > ${TAXON}_CDS.blank

#Now we join this file to our original .rc file and the result includes entries in the blank file that were not present in the original .rc

join -t$'\t' -a2 <(awk '{OFS="\t"; print $1"!"$2, $0}' ${TAXON}.rc | sort -k1,1) <(awk '{OFS="\t"; print $1"!"$2,$0}' ${TAXON}_CDS.blank | sort -k1,1) | cut --complement -f1 |sort -k1,1 -k2,2n > ${TAXON}.rc2

#I think this last part worked, but it looks like it contains both the .rc and blank file. I'll proceed and see what happens.

#Arrange info in .rc2 into a .tsv file
awk '{split($6,a,":");split($7,c,":");split($8,g,":");split($9,t,":");b=$2%3;printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, b, $3, a[2]+c[2]+g[2]+t[2], a[2], c[2], g[2], t[2]}' ${TAXON}.rc2 > ${TAXON}.tsv

#Check average read depth and exclude genes with an average depth of coverage <20 reads.

