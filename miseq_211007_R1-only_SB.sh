#!/bin/bash

################
# set up
################

cd /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-opt-new

# list files in folder
# extract part of name before second underscore and find unique
ls input | cut -d_ -f1,2 | sort | uniq

# store list of names in bash variable
samples=$(ls input | cut -d_ -f1,2 | sort | uniq)
# check contents
echo $samples

# count number of sequences across all files in folder
cd input
gzcat *.fastq.gz | grep -c "^@M00"  
# 3,084,014 (BW)
# 2,197,228 (SB - no b2 fish)
# 2,265,595 (SB - b2 fish, still less than BW though - why? NTCs?)
# 2,321,597 (B2 + NTCs)

for s in $samples;
do
echo "${s}_L001_R1_001.fastq.gz" \
count=$(gunzip -c ${s}_L001_R1_001.fastq.gz | grep -c "^@M00" ) \
>> SeqNumber_Smal-corII_211007_input.txt;
done


cd ..
mkdir output


################
# Primer removal 
################
# do not anchor (^) primers
# because R1 only and 
# amplicon length ~ 90bp
# Smal-corII: GGCGCTTGTAACGCCAAGGTAGTGGGTTCGATCCCCGGGACCACCCATACACAAAAATGTATGCACGCATGACTGTAAGTCGCTTTGGAT
# use linked (Fwd...rcRev) 
# expected amplicon length 90bp (Hamada et al. 1997)
# set minimum length (-m) to 70
#can play around with this 

for s in $samples;
do
cutadapt -g TAGCTCAGCTGGTAGAGCAC...AAAAGCGTCTGCTAAATGGCA  \
-o output/${s}_L001_R1.fastq.gz --discard-untrimmed -m 70 \
input/${s}_L001_R1_001.fastq.gz;
done


cd ../output
gzcat *.fastq.gz | grep -c "^@M00"  
#1,311,257 (BW)
#1,219,463 (SB - no b2 fish)
#1,284,790 (SB - b2 fish, still less than BW though - why?)
#1,292,472 (B2 + NTCs)

# count number of sequences per file and output as list
# must be in the folder where you want to count sequences

for s in $samples;
do
echo "${s}_L001_R1_001.fastq.gz" \
count=$(gunzip -c ${s}_L001_R1.fastq.gz | grep -c "^@M00" ) \
>> SeqNumber_Smal-corII_211007_primersgone.txt;
done

#####################
# Sequence Quality
#####################
cd ../output

# IF fusing rather than merging do this after quality filtering
# make sure in folder 'trimmed'
mkdir qc

# Check sequence quality
# Run fastqc on each file in input
ls *.fastq.gz | parallel 'fastqc {}'

#Move qc outputs
mv /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-opt-new/output/*.html /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-opt-new/output/qc
mv /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-opt-new/output/*.zip /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-opt-new/output/qc

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

#can tighten up these length paramters based on what I'm seeing with the gDNA

for s in $samples;
do
cutadapt -o trimmed/${s}_R1_trimmed.fastq.gz -q 26 -m 70 -M 105 \
./${s}_L001_R1.fastq.gz;
done

#look at the output to make sure the adaptor is gone

#convert fastq to fasta
mkdir fasta

for s in $samples;
do
vsearch --fastx_filter trimmed/${s}_R1_trimmed.fastq.gz --fastaout fasta/${s}_R1_trimmed.fasta
done

################
# Dereplication
################

mkdir derep
# derep individuals samples

for s in $samples;
do
vsearch --derep_fulllength fasta/${s}_R1_trimmed.fasta --sizeout --minuniquesize 10 --output derep/derep_${s}_R1_trimmed.fasta;
done 











## trying it the old way but with these parameters to see if it's the same ##

################ Input data ################

# list files in folder
# extract part of name before second underscore and find unique
ls input | cut -d_ -f1,2 | sort | uniq

# store list of names in bash variable and check contents
samples=$(ls input | cut -d_ -f1,2 | sort | uniq)
echo $samples


# count number of sequences across all files in folder
cd input
gzcat *.fastq.gz | grep -c "^@M00" 
#2,197,228

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

# ~20K-70K reads/sample



################ Primer removal ################

cd .. 

mkdir output1

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
-o output1/${s}_L001_R1.fastq.gz --discard-untrimmed -m 70 \
input/${s}_L001_R1_001.fastq.gz;
done


# count number of sequences across all files in folder
cd input
gzcat *.fastq.gz | grep -c "^@M00"  
# 3,084,014 (BW)
# 2,197,228 (SB - no b2 fish)
# 2,265,595 (SB - b2 fish, still less than BW though - why? NTCs?)
# 2,321,597 (B2 + NTCs)

cd ../output1
gzcat *.fastq.gz | grep -c "^@M00"  
#1,311,257 (BW)
#1,219,463 (SB - no b2 fish)
#1,284,790 (SB - b2 fish, still less than BW though - why?)
#1,292,472 (B2 + NTCs)

# count number of sequences per file and output as list
# must be in the folder where you want to count sequences
cd ../input

for s in $samples;
do
echo "${s}_L001_R1_001.fastq.gz" \
count=$(gunzip -c ${s}_L001_R1_001.fastq.gz | grep -c "^@M00" ) \
>> SeqNumber_Smal-corII_211007_input.txt;
done



#####################
# Sequence Quality
#####################
cd ../output1

# IF fusing rather than merging do this after quality filtering
# make sure in folder 'trimmed'
mkdir qc

# Check sequence quality
# Run fastqc on each file in input
ls *.fastq.gz | parallel 'fastqc {}'

#Move qc outputs
mv /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-opt-new/output1/*.html /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-opt-new/output1/qc
mv /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-opt-new/output1/*.zip /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-opt-new/output1/qc

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
#1,274,208



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
#1,274,208




################ convert to fasta ################
# make sure in ASV folder

vsearch --fastx_filter trimmed_concatenated.fastq --fastaout trimmed_concatenated.fasta


#count seqs - this way is not informative as to the depth/fish, only total
grep -c ">M00" trimmed_concatenated.fasta
#1,274,207 (one less bc of the .txt file not getting written over)



################ Dereplication ################

# header of each unique sequence records number of copies in dataset

vsearch --derep_fulllength trimmed_concatenated.fasta --output derep.fasta --sizeout --minuniquesize 10 --relabel uniq

#count seqs
vsearch --search_exact trimmed_concatenated.fasta -db derep.fasta -otutabout derep.tsv

#these two methods give me the same results :) 
#will use the one with concatenation so can keep using my downstream scripts

