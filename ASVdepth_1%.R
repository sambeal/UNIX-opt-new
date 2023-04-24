#analysis of the ASVs returned from the UNIX pipeline when run with default parameters except:
#sequence length restricted to 70-105bp (as determined from multiqc report)
#phred score 26 enforced during trimming
#sequences only taken through the dereplication stage

#Parse through true ASVs to find only those with depths >=10% per fish

#set working directory----
getwd()

#"/Users/samanthabeal"
setwd("Documents/MSc/Bioinformatics")


#read in data----
true <- read.csv("UNIX-opt-new/TrueASVs.csv", header = TRUE)

#count # unique ASVs
unique <- as.data.frame(unique(true$ASV))
#8891

#add total depth
true$Total.depth <- apply(true[,2:4],1,sum)

#rearrange
library(dplyr)
true <- true %>% relocate(Total.depth, .before = Fish)


#add in seqs----
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
true_seqs <- true %>% 
  inner_join(UNIXseqs, by = c("ASV" = "ASV"))

#rearrange
true_seqs <- true_seqs %>% relocate(Sequence, .after = ASV)

#add seq lengths
true_seqs$Sequence_length <- str_length(true_seqs$Sequence)
true_seqs <- true_seqs %>% relocate(Sequence_length, .after = Sequence)

#pull out only ASVs present with >=1% depth per fish----
#10% seems too high --> A1190222 @ 10% = 1 ASV

#190222-A1----
A1190222 <- subset(true_seqs, Fish == "A1190222", c("ASV", "Sequence","Sequence_length","Total.depth","Fish"))
sum(A1190222$Total.depth)
#39,589

#find cut off of this fish
A1190222_1perc <- subset(A1190222, Total.depth >= (sum(A1190222$Total.depth)*0.01))
sum(A1190222_1perc$Total.depth)
#12,122

#190222-A2----
A2190222 <- subset(true_seqs, Fish == "A2190222", c("ASV", "Sequence","Sequence_length","Total.depth","Fish"))
sum(A2190222$Total.depth)
#36,688

#find cut off of this fish
A2190222_1perc <- subset(A2190222, Total.depth >= (sum(A2190222$Total.depth)*0.01))
sum(A2190222_1perc$Total.depth)
#12,967

#190222-H1----
H1190222 <- subset(true_seqs, Fish == "H1190222", c("ASV", "Sequence","Sequence_length","Total.depth","Fish"))
sum(H1190222$Total.depth)
#807

#find cut off of this fish
H1190222_1perc <- subset(H1190222, Total.depth >= (sum(H1190222$Total.depth)*0.01))
sum(H1190222_1perc$Total.depth)
#789

#190222-M1----
M1190222 <- subset(true_seqs, Fish == "M1190222", c("ASV", "Sequence","Sequence_length","Total.depth","Fish"))
sum(M1190222$Total.depth)
#38,041

#find cut off of this fish
M1190222_1perc <- subset(M1190222, Total.depth >= (sum(M1190222$Total.depth)*0.01))
sum(M1190222_1perc$Total.depth)
#13,990

#190918-A1----
A1190918 <- subset(true_seqs, Fish == "A1190918", c("ASV", "Sequence","Sequence_length","Total.depth","Fish"))
sum(A1190918$Total.depth)
#13,110

#find cut off of this fish
A1190918_1perc <- subset(A1190918, Total.depth >= (sum(A1190918$Total.depth)*0.01))
sum(A1190918_1perc$Total.depth)
#11,960

write_csv(A1190918, "UNIX-opt-new/A1190918_true.csv")
write_csv(A1190918_1perc, "UNIX-opt-new/A1190918_1perc.csv")

#190918-A3----
A3190918 <- subset(true_seqs, Fish == "A3190918", c("ASV", "Sequence","Sequence_length","Total.depth","Fish"))
sum(A3190918$Total.depth)
#40,209

#find cut off of this fish
A3190918_1perc <- subset(A3190918, Total.depth >= (sum(A3190918$Total.depth)*0.01))
sum(A3190918_1perc$Total.depth)
#12,815

#190918-A4----
A4190918 <- subset(true_seqs, Fish == "A4190918", c("ASV", "Sequence","Sequence_length","Total.depth","Fish"))
sum(A4190918$Total.depth)
#37,951

#find cut off of this fish
A4190918_1perc <- subset(A4190918, Total.depth >= (sum(A4190918$Total.depth)*0.01))
sum(A4190918_1perc$Total.depth)
#12,545

#190918-C4----
C4190918 <- subset(true_seqs, Fish == "C4190918", c("ASV", "Sequence","Sequence_length","Total.depth","Fish"))
sum(C4190918$Total.depth)
#45,855

#find cut off of this fish
C4190918_1perc <- subset(C4190918, Total.depth >= (sum(C4190918$Total.depth)*0.01))
sum(C4190918_1perc$Total.depth)
#14,636

#190918-D4----
D4190918 <- subset(true_seqs, Fish == "D4190918", c("ASV", "Sequence","Sequence_length","Total.depth","Fish"))
sum(D4190918$Total.depth)
#42,193

#find cut off of this fish
D4190918_1perc <- subset(D4190918, Total.depth >= (sum(D4190918$Total.depth)*0.01))
sum(D4190918_1perc$Total.depth)
#13,127

#190918-E7----
E7190918 <- subset(true_seqs, Fish == "E7190918", c("ASV", "Sequence","Sequence_length","Total.depth","Fish"))
sum(E7190918$Total.depth)
#37,831

#find cut off of this fish
E7190918_1perc <- subset(E7190918, Total.depth >= (sum(E7190918$Total.depth)*0.01))
sum(E7190918_1perc$Total.depth)
#12,423

#190918-F2----
F2190918 <- subset(true_seqs, Fish == "F2190918", c("ASV", "Sequence","Sequence_length","Total.depth","Fish"))
sum(F2190918$Total.depth)
#45,544

#find cut off of this fish
F2190918_1perc <- subset(F2190918, Total.depth >= (sum(F2190918$Total.depth)*0.01))
sum(F2190918_1perc$Total.depth)
#14,931

#190922-B2----
B2190922 <- subset(true_seqs, Fish == "B2190922", c("ASV", "Sequence","Sequence_length","Total.depth","Fish"))
sum(B2190922$Total.depth)
#739

#find cut off of this fish
B2190922_1perc <- subset(B2190922, Total.depth >= (sum(B2190922$Total.depth)*0.01))
sum(B2190922_1perc$Total.depth)
#728

#190922-D2----
D2190922 <- subset(true_seqs, Fish == "D2190922", c("ASV", "Sequence","Sequence_length","Total.depth","Fish"))
sum(D2190922$Total.depth)
#31,545

#find cut off of this fish
D2190922_1perc <- subset(D2190922, Total.depth >= (sum(D2190922$Total.depth)*0.01))
sum(D2190922_1perc$Total.depth)
#10,811

#190922-F2----
F2190922 <- subset(true_seqs, Fish == "F2190922", c("ASV", "Sequence","Sequence_length","Total.depth","Fish"))
sum(F2190922$Total.depth)
#37,170

#find cut off of this fish
F2190922_1perc <- subset(F2190922, Total.depth >= (sum(F2190922$Total.depth)*0.01))
sum(F2190922_1perc$Total.depth)
#12,931

#190922-G2----
G2190922 <- subset(true_seqs, Fish == "G2190922", c("ASV", "Sequence","Sequence_length","Total.depth","Fish"))
sum(G2190922$Total.depth)
#39,750

#find cut off of this fish
G2190922_1perc <- subset(G2190922, Total.depth >= (sum(G2190922$Total.depth)*0.01))
sum(G2190922_1perc$Total.depth)
#13,326

#190922-H2----
H2190922 <- subset(true_seqs, Fish == "H2190922", c("ASV", "Sequence","Sequence_length","Total.depth","Fish"))
sum(H2190922$Total.depth)
#41,629

#find cut off of this fish
H2190922_1perc <- subset(H2190922, Total.depth >= (sum(H2190922$Total.depth)*0.01))
sum(H2190922_1perc$Total.depth)
#14,025

#190922-J2----
J2190922 <- subset(true_seqs, Fish == "J2190922", c("ASV", "Sequence","Sequence_length","Total.depth","Fish"))
sum(J2190922$Total.depth)
#40,693

#find cut off of this fish
J2190922_1perc <- subset(J2190922, Total.depth >= (sum(J2190922$Total.depth)*0.01))
sum(J2190922_1perc$Total.depth)
#12,781

#NTC----
NTC <- subset(true_seqs, Fish == "NTC", c("ASV", "Sequence","Sequence_length","Total.depth","Fish"))
sum(NTC$Total.depth)
#0

#find cut off of this fish
NTC_1perc <- subset(NTC, Total.depth >= (sum(NTC$Total.depth)*0.01))
sum(NTC_1perc$Total.depth)
#0


#put it all together----
all_1perc <- rbind(A1190222_1perc,A2190222_1perc, H1190222_1perc,M1190222_1perc,
                   A1190918_1perc,A3190918_1perc,A4190918_1perc,C4190918_1perc,
                   D4190918_1perc, E7190918_1perc,F2190918_1perc,B2190922_1perc,
                   D2190922_1perc, F2190922_1perc,G2190922_1perc,H2190922_1perc,
                   J2190922_1perc,NTC_1perc)




#make file of unique ASV name, seq, and total depth across all fish
ASVs <- aggregate(Total.depth ~ ASV + Sequence, data=all_1perc, FUN = sum)


#add sequence length back in
ASVs$Sequence_length <- str_length(ASVs$Sequence)
ASVs <- ASVs %>% relocate(Sequence_length, .after = Sequence)


#want to add a column for the number of ASVs per fish
# OR 
# the number of fish each ASV was in
ASVsperFish <- aggregate(ASV ~ Fish , data=all_1perc, FUN = length)
FishperASV <- aggregate(Fish ~ ASV , data=all_1perc, FUN = length)


#add FishperASV to ASV file by "ASV"
library(tidyverse)
ASVs_fish <- ASVs %>% 
  inner_join(FishperASV , by = c("ASV" = "ASV"))
names(ASVs_fish)[4] <- "Num.Fish"
ASVs_fish$Sequence_length <- str_length(ASVs_fish$Sequence)


#save output----
write_csv(all_1perc,"UNIX-opt-new/ASVdepth_1percent.csv")



#make fasta for downstream analysis----
#make data frame of only ASV and seq then convert to .fasta
UNIXfasta <- ASVs[c(1:2)]
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
            file = "UNIX-opt-new/ASVdepth_1percent.fasta", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)




#alig the sequences----

# load the sequences from the file
# change "DNA" to "RNA" or "AA" if necessary
library("Biostrings")
library("here")
#"here() starts at /Users/samanthabeal/Documents/MSc/Bioinformatics"

ASVs1perc <- readDNAStringSet(here("UNIX-opt-new/ASVdepth_1percent.fasta"), format = "fasta")

# look at some of the sequences (optional)
ASVs1perc

# nucleotide sequences need to be in the same orientation
# if they are not, then they can be reoriented (optional)
library("DECIPHER")
ASVs1perc <- OrientNucleotides(ASVs1perc)

# perform the alignment
alignedASVs <- AlignSeqs(ASVs1perc)

# view the alignment in a browser (optional)
BrowseSeqs(alignedASVs, highlight=0)

# write the alignment to a new FASTA file
writeXStringSet(alignedASVs,
                file="UNIX-opt-new/alignedASVs1perc.fasta")

#plot----
#histogram----
hist(ASVs_fish$Sequence_length, main="Distribution of SINE ASVs from gDNA",
     xlab="Length of ASV seqeuence")

hist(ASVs_fish$Num.Fish, main="Number of fish containing the ASV",
     xlab="Number of fish per ASV")

#(hamming) distance matrix----
library(matrixcalc)
ASVs1_diffs <- as.matrix(stringDist(alignedASVs,method="hamming"))

# Visualize distance matrix
library(ggplot2)

#Concert matrix to long form to work with ggplot heatmap
ASVs1_diffs_long <- ASVs1_diffs %>% as.data.frame() %>%
  rownames_to_column("ASV") %>%
  pivot_longer(-c(ASV), names_to="ASV2",values_to="distance")

ggplot(ASVs1_diffs_long, aes(x=ASV,y=ASV2,fill=distance))+
  geom_raster()+
  scale_fill_gradient(low="white",high="darkred")+
  theme(axis.text.x=element_text(angle=90))+
  ggtitle("Base Differences Among UNIX SINE ASVs")


#analyze which fish per ASV and which ASVs per fish from ASVdepth_1percent.csv----


#make haplotype network----
library("ape")
library("pegas")

UNIXdata <-read.dna("UNIX-opt-new/alignedASVs1perc.fasta", format="fasta")
UNIXhaplo <- haplotype(UNIXdata)
UNIXhaplo

# Number of haplotypes:25, Sequence length: 93
UNIXhaplodist <- dist.dna(UNIXhaplo, "N")
UNIXnetwork <- rmst(UNIXhaplodist, quiet = TRUE)
UNIXnetwork

#Haplotype network with:
#25 haplotypes
#90 links
#link lengths between 0 and 2 steps
plot(UNIXnetwork)
plot(UNIXnetwork, fast = TRUE)
#dont like how this network looks

#need to convert the aligned.fasta to .nex for using POPart to make the tree
#did a command line conversion
#seqmagick convert --output-format nexus --alphabet dna Documents/MSc/Bioinformatics/UNIX-opt-new/alignedASVs1perc.fasta Documents/MSc/Bioinformatics/UNIX-opt-new/alignedASVs1perc.nex








