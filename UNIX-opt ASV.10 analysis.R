#analysis of the ASVs returned from the UNIX pipeline when run with default parameters except:
#sequence length restricted to 70-105bp (as determined from multiqc report)
#phred score 26 enforced during trimming
#sequences only taken through the dereplication stage

#Parse through true ASVs to find only those >10 occurances----

#set working directory
getwd()
#"/Users/samanthabeal"
setwd("Documents/MSc/Bioinformatics")


#read in data
true <- read.csv("UNIX-opt-new/TrueASVs.csv", header = TRUE)

#count # unique ASVs
unique <- as.data.frame(unique(true$ASV))
#8891

#count unique ASVs per Fish
aggregate(ASV ~ Fish, data=true, FUN=function(ASV) length(unique(ASV)))
#      Fish  ASV
#1  A1190222 2448
#2  A1190918   89
#3  A2190222 1996
#4  A3190918 2364
#5  A4190918 2325
#6  C4190918 2722
#7  D2190922 1823
#8  D4190918 2579
#9  E7190918 2275
#10 F2190918 2651
#11 F2190922 2050
#12 G2190922 2419
#13 H1190222   15
#14 H2190922 2220
#15 J2190922 2441
#16 M1190222 1881


#pull out only ASVs occuring >10x in at least one PCR rep
true10 <- subset(true, PCR.1 >= 10 | PCR.2 >= 10 | PCR.3 >= 10)
#3406

#count unique ASVs >=10 depth in at least PCR rep per Fish
aggregate(ASV ~ Fish, data=true10, FUN=function(ASV) length(unique(ASV)))
#Fish ASV
#1  A1190222 242
#2  A1190918  40
#3  A2190222 237
#4  A3190918 244
#5  A4190918 225
#6  C4190918 259
#7  D2190922 192
#8  D4190918 249
#9  E7190918 232
#10 F2190918 267
#11 F2190922 234
#12 G2190922 215
#13 H1190222   8
#14 H2190922 257
#15 J2190922 233
#16 M1190222 272



#add total depth
true10$Total.depth <- apply(true10[,2:4],1,sum)

#rearrange
library(dplyr)
true10 <- true10 %>% relocate(Total.depth, .before = Fish)

#add in seqs
library (devtools)
source_url("https://raw.githubusercontent.com/lrjoshi/FastaTabular/master/fasta_and_tabular.R")

FastaToTabular("UNIX-opt-new/output1/ASV/derep.fasta")
#The output will be stored as dna_table.csv in the current directory ("Bioinformatics")

#read in data
UNIXseqs <- read.csv("dna_table.csv")

#reformat
UNIXseqs <- UNIXseqs[-1]
names(UNIXseqs)[1] <- "ASV"
names(UNIXseqs)[2] <- "Sequence"

#";.*" says: replace semicolon (;) and every character after that (.*), with nothing "".
UNIXseqs$ASV <- gsub(";.*", "", UNIXseqs$ASV)
#">" says: replace arrow (>) with nothing "".
UNIXseqs$ASV <- gsub(">", "", UNIXseqs$ASV)

#add data frames together
library(tidyverse)
true10_seqs <- true10 %>% 
  inner_join(UNIXseqs, by = c("ASV" = "ASV"))

#rearrange
true10_seqs <- true10_seqs %>% relocate(Sequence, .after = ASV)


#pull out only unique ASVs
true10unique <- as.data.frame(unique(true10_seqs$ASV))
names(true10unique)[1] <- "ASV"

#make file of ASV name, seq, and depth
ASVs <- true10_seqs[c(1,2,6)]

true_seqs_full <- true10unique %>% 
  inner_join(ASVs, by = c("ASV" = "ASV"))

true_seqs <- aggregate(Total.depth ~ ASV + Sequence, data = true_seqs_full, FUN = sum)

#add seq lengths
true_seqs$Sequence_length <- str_length(true_seqs$Sequence)


#save output
write_csv(true_seqs,"UNIX-opt-new/ASVdepth10.csv")

hist(true_seqs$Sequence_length, main="Distribution of SINE ASVs from gDNA",
     xlab="Length of ASV seqeuence")


#make fasta for downstream analysis----
#make data frame of only ASV and seq then convert to .fasta
UNIXfasta <- true_seqs[c(1:2)]
names(UNIXfasta)[1] <- "ASV"
UNIXfasta$ASV <- paste0(">", UNIXfasta$ASV)

library(tidyverse)

UNIXall_print <-
  UNIXfasta %>%
  select(ASV,Sequence) %>%
  rowwise() %>%
  pivot_longer(ASV:Sequence) %>%
  select(-name)

#allseqs
write.table(UNIXall_print, 
            file = "UNIX-opt-new/ASVdepth10.fasta", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)




#alig the sequences----

# load the sequences from the file
# change "DNA" to "RNA" or "AA" if necessary
library("Biostrings")
library("here")
#"here() starts at /Users/samanthabeal/Documents/MSc/Bioinformatics"

ASVs10 <- readDNAStringSet(here("UNIX-opt-new/ASVdepth10.fasta"), format = "fasta")

# look at some of the sequences (optional)
ASVs10

# nucleotide sequences need to be in the same orientation
# if they are not, then they can be reoriented (optional)
library("DECIPHER")
ASVs10 <- OrientNucleotides(ASVs10)

# perform the alignment
alignedASVs <- AlignSeqs(ASVs10)

# view the alignment in a browser (optional)
BrowseSeqs(alignedASVs, highlight=0)

# write the alignment to a new FASTA file
writeXStringSet(alignedASVs,
                file="UNIX-opt-new/alignedASVdepth10.fasta")







# Find differences among sequences (Hamming distance)----
library(matrixcalc)
ASVs10_diffs <- as.matrix(stringDist(alignedASVs,method="hamming"))

# Visualize distance matrix
library(ggplot2)

#Concert matrix to long form to work with ggplot heatmap
ASVs10_diffs_long <- ASVs10_diffs %>% as.data.frame() %>%
  rownames_to_column("ASV") %>%
  pivot_longer(-c(ASV), names_to="ASV2",values_to="distance")

ggplot(ASVs10_diffs_long, aes(x=ASV,y=ASV2,fill=distance))+
  geom_raster()+
  scale_fill_gradient(low="white",high="darkred")+
  theme(axis.text.x=element_text(angle=90))+
  ggtitle("Base Differences Among UNIX SINE ASVs")

#make haplotype network----
library("ape")
library("pegas")

UNIXdata <-read.dna("UNIX-opt-new/alignedASVdepth10.fasta", format="fasta")
UNIXhaplo <- haplotype(UNIXdata)
UNIXhaplo

# Number of haplotypes:405, Sequence length: 101 

UNIXhaplodist <- dist.dna(UNIXhaplo, "N")
UNIXnetwork <- rmst(UNIXhaplodist, quiet = TRUE)
UNIXnetwork

#Haplotype network with:
#405 haplotypes
#22028 links
#link lengths between 0 and 2 steps
plot(UNIXnetwork)
plot(UNIXnetwork, fast = TRUE)







