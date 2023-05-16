#!/bin/bash

################ Set up ################

cd /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-opt-new/eDNA/release_02


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
mkdir output

# do not anchor (^) primers
# amplicon length ~ 90bp
# Smal-corII: GGCGCTTGTAACGCCAAGGTAGTGGGTTCGATCCCCGGGACCACCCATACACAAAAATGTATGCACGCATGACTGTAAGTCGCTTTGGAT
# expected amplicon length 90bp (Hamada et al. 1997)
# set minimum length (-m) to 70


for s in $samples;
do
	cutadapt -g TAGCTCAGCTGGTAGAGCAC -G TGCCATTTAGCAGACGCTTTT \
	-o output/${s}_L001_R1.fastq.gz -p output/${s}_L001_R2.fastq.gz --discard-untrimmed -m 70 \
	input/${s}_L001_R1_001.fastq.gz input/${s}_L001_R2_001.fastq.gz;
done


cd output
gzcat *.fastq.gz | grep -c "^@M00"



##################### Sequence Quality #####################

# IF fusing rather than merging do this after quality filtering
# make sure in folder 'output'
mkdir qc

# Check sequence quality
# Run fastqc on each file in input
ls *.fastq.gz | parallel 'fastqc {}'

#Move qc outputs
mv /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-opt-new/eDNA/release_02/output/*.html /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-opt-new/eDNA/release_02/output/qc
mv /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-opt-new/eDNA/release_02/output/*.zip /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-opt-new/eDNA/release_02/output/qc

cd qc
multiqc .

#shingle qc = very poor. 31/54 samples failed seq quality histogram
#release_02 qc = very poor. 32/50 samples failed seq quality histogram


################ Merge pairs (iF forward and reverse sequences overlap) ################

# Smal cor-II PCR product length: 131
# 131 - (18 + 18) = 95 -> OVERLAP!! try setting -v 60
# -q sets minimum quality
#(We generally set the quality threshold (-q) for trimming low quality ends to 26 - learn to metabarcode)
# -v sets minimum overlap
#(minimum overlap (-v) to a value around 20-30 bases less than the presumed overlap length)

cd .. 
mkdir merged
for s in $samples; do pear -f ${s}_L001_R1.fastq.gz -r ${s}_L001_R2.fastq.gz -o merged/$s -q 26 -v 18; done
#default min overlap is 10


# move unassembled and discarded.fastq in new folder
cd merged
mkdir unassembled_discarded
mv /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-opt-new/eDNA/release_02/output/merged/*discarded.fastq /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-opt-new/eDNA/release_02/output/merged/unassembled_discarded
mv /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-opt-new/eDNA/release_02/output/merged/*unassembled*.fastq /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-opt-new/eDNA/release_02/output/merged/unassembled_discarded

# remove assembled from file name
rename -e "s/assembled\.//" *


 
################ Quality Trim Sequences ################

cd ../ #(back to output)
mkdir trimmed

# quality (-q)
# minimum length (-m)
# trim based on qc report
# maximum length (-M)

cd merged
for s in $samples;
do
cutadapt -o ../trimmed/${s}_trimmed.fastq -q 26 -m 70 -M 105 \
./${s}.fastq;
done




################ Concatenate Samples ################

cd ..
mkdir ASV #(in output)

cd trimmed


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




################ Convert to fasta ################

# make sure in ASV folder
vsearch --fastx_filter trimmed_concatenated.fastq --fastaout trimmed_concatenated.fasta

#count seqs - make sure didn't lose any during the conversion
grep -c ">M00" trimmed_concatenated.fasta




################ Dereplication ################

# header of each unique sequence records number of copies in dataset

vsearch --derep_fulllength trimmed_concatenated.fasta --output derep.fasta --sizeout --minuniquesize 10 --relabel uniq

#count seqs
vsearch --search_exact trimmed_concatenated.fasta -db derep.fasta -otutabout derep.tsv




#unique ASVs are lost when taking samples through the rest of the pipeline. Stop here




################ Denoising ################

vsearch --cluster_unoise derep.fasta --minsize 8 --unoise_alpha 2 --centroids denoised.fasta



################ Chimera Filtering ################
 
vsearch --uchime3_denovo denoised.fasta --nonchimeras nochim.fasta




####################### Mapping Reads to ASVs #######################

vsearch --search_exact trimmed_concatenated.fasta -db nochim.fasta -otutabout nochim.tsv





