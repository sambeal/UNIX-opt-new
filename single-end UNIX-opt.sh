#!/bin/bash

################
# set up
################

cd /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-opt-new


################ Input data ################
cd gDNA

# list files in folder
# extract part of name before second underscore and find unique
ls input | cut -d_ -f1,2 | sort | uniq

# store list of names in bash variable and check contents
samples=$(ls input | cut -d_ -f1,2 | sort | uniq)
echo $samples


# count number of sequences across all files in folder
cd input
gzcat *.fastq.gz | grep -c "^@M00" 


# count number of sequences per file and output as list
# must be in the folder where you want to count sequences
#initial input files
#gunzip = g unzip (opposit of this is gzip)

for s in $samples;
do
echo "${s}_L001_R1_001.fastq.gz" \
count=$(gunzip -c ${s}_L001_R1_001.fastq.gz | grep -c "^@M00" ) \
>> SeqNumber_input.txt;
done



################ Primer removal ################

cd .. 

mkdir output

# do not anchor (^) primers
# because R1 only and 
# amplicon length ~ 90bp
# Smal-corII: GGCGCTTGTAACGCCAAGGTAGTGGGTTCGATCCCCGGGACCACCCATACACAAAAATGTATGCACGCATGACTGTAAGTCGCTTTGGAT
# use linked (Fwd...rcRev) 
# expected amplicon length 90bp (Hamada et al. 1997)
# set minimum length (-m) to 70


for s in $samples;
do
cutadapt -g TAGCTCAGCTGGTAGAGCAC...AAAAGCGTCTGCTAAATGGCA  \
-o output/${s}_L001_R1.fastq.gz --discard-untrimmed -m 70 \
input/${s}_L001_R1_001.fastq.gz;
done


# count number of sequences across all files in folder
cd input
gzcat *.fastq.gz | grep -c "^@M00"  
#2321597

cd ../output
gzcat *.fastq.gz | grep -c "^@M00"  
#1292472

for s in $samples;
do
echo "${s}_L001_R1_001.fastq.gz" \
count=$(gunzip -c ${s}_L001_R1.fastq.gz | grep -c "^@M00" ) \
>> SeqNumber_primerremoval.txt;
done



#####################
# Sequence Quality
#####################
#make sure in output

# IF fusing rather than merging do this after quality filtering
# make sure in folder 'trimmed'
mkdir qc

# Check sequence quality
# Run fastqc on each file in input
ls *.fastq.gz | parallel 'fastqc {}'

#Move qc outputs
mv /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-opt-new/gDNA/output/*.html /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-opt-new/gDNA/output/qc
mv /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-opt-new/gDNA/output/*.zip /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-opt-new/gDNA/output/qc

cd qc
multiqc .

# quality drops at ~ 105bp
# very messy quality histogram


#########################
# Quality Trim Sequences
########################
cd ../
mkdir trimmed

# quality (-q)
# minimum length (-m)
# trim based on qc report
# maximum length (-M)

for s in $samples;
do
cutadapt -o trimmed/${s}_R1_trimmed.fastq.gz -q 26 -m 70 -M 105 \
./${s}_L001_R1.fastq.gz;
done

cd trimmed
for s in $samples;
do
            echo "${s}_L001_R1_001.fastq.gz" \
            count=$(gunzip -c ${s}_R1_trimmed.fastq.gz | grep -c "^@M00" ) \
       >> SeqNumber_trimmed.txt;
done

gzcat *.fastq.gz | grep -c "^@M00" 
#1274208



################ Concatenate Samples ################
cd ..
mkdir ASV #(in output)

cd trimmed
gunzip *.fastq.gz

# add sample name to sequences
# combine all sequences in one file
for f in *;
do                         
        sed -e "s/\(^@M00.*\) .*$/\1;sample=${f%.*};/" $f \
        >> ../ASV/trimmed_concatenated.fastq;
done

#count seqs
cd ../ASV
grep -c "M00" trimmed_concatenated.fastq
#1274208




################ convert to fasta ################
# make sure in ASV folder

vsearch --fastx_filter trimmed_concatenated.fastq --fastaout trimmed_concatenated.fasta


#count seqs - this way is not informative as to the depth/fish, only total
grep -c ">M00" trimmed_concatenated.fasta
#1274207 (one less bc cannot count seq count file)



################ Dereplication ################

# header of each unique sequence records number of copies in dataset

vsearch --derep_fulllength trimmed_concatenated.fasta --output derep.fasta --sizeout --minuniquesize 10 --relabel uniq

#count seqs
vsearch --search_exact trimmed_concatenated.fasta -db derep.fasta -otutabout derep.tsv




