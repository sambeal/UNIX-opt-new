#analysis of eDNA ASVs
#230509 = Aquatron
#230510 = Shingle (with ASP pos con and without)
#230512 = AWR000_02 (with ASP pos con removed)

#load data
all_1perc_seqs <- read.csv("eDNA/aquatron/230509/1%_ASVs.csv")

#assess the retained ASVs
all_1perc_seqs$ASV <- as.factor(all_1perc_seqs$ASV)

#plotting----
#fish/ASV tile matrix----
#pull out only sample and ASV
sample_ASVs <- all_1perc_seqs[c(1,7)]
sample_ASVs$ASV <- as.factor(sample_ASVs$ASV)
sample_ASVs$Sample <- as.factor(sample_ASVs$Sample)

uniqueASV <- as.data.frame(unique(sample_ASVs$ASV))
names(uniqueASV)[1] <- "ASV"
#8 unique ASVs (Aquatron)
#73 (Shingle)

#make dataframe of unique fish for the matrix plot
uniqueSample <- as.data.frame(unique(sample_ASVs$Sample))
names(uniqueSample)[1] <- "Sample"

#subset into a dataframe for each ASV (less ASVs than fish)
ASV_split <- split(all_1perc_seqs, all_1perc_seqs$ASV)
str(ASV_split)
uniq1 <- ASV_split[[1]]
uniq2 <- ASV_split[[2]]
uniq3 <- ASV_split[[3]]
uniq4 <- ASV_split[[4]]
uniq5 <- ASV_split[[5]]
uniq6 <- ASV_split[[6]]
uniq7 <- ASV_split[[7]]
uniq8 <- ASV_split[[8]]

#check if fish in in each ASV df
#uniq1
uniq1_sample <- as.data.frame(uniqueSample$Sample %in% uniq1$Sample)
names(uniq1_sample)[1] <- "uniq1"

#uniq2
uniq2_sample <- as.data.frame(uniqueSample$Sample %in% uniq2$Sample)
names(uniq2_sample)[1] <- "uniq2"

#uniq3
uniq3_sample <- as.data.frame(uniqueSample$Sample %in% uniq3$Sample)
names(uniq3_sample)[1] <- "uniq3"

#uniq4
uniq4_sample <- as.data.frame(uniqueSample$Sample %in% uniq4$Sample)
names(uniq4_sample)[1] <- "uniq4"

#uniq5
uniq5_sample <- as.data.frame(uniqueSample$Sample %in% uniq5$Sample)
names(uniq5_sample)[1] <- "uniq5"

#uniq6
uniq6_sample <- as.data.frame(uniqueSample$Sample %in% uniq6$Sample)
names(uniq6_sample)[1] <- "uniq6"

#uniq7
uniq7_sample <- as.data.frame(uniqueSample$Sample %in% uniq7$Sample)
names(uniq7_sample)[1] <- "uniq7"

#uniq8
uniq8_sample <- as.data.frame(uniqueSample$Sample %in% uniq8$Sample)
names(uniq8_sample)[1] <- "uniq8"

#add together
ASVperSample <- cbind(uniq1_sample,uniq2_sample,uniq3_sample,
                    uniq4_sample, uniq5_sample, uniq6_sample,
                    uniq7_sample, uniq8_sample)

#add sample names back in
rownames(ASVperSample) <- uniqueSample[,1]



#save
library(readr)
write_csv(ASVperSample, "eDNA/aquatron/230509/true_false_matrix.csv")

#change from true/false to 1/0 by making a matrix where true=1, false=0
matrix <- data.matrix(ASVperSample)

#plot
library(ggplot2)
library(reshape2)

#Convert wide to long
matrix.long <- melt(matrix)
matrix.long$value <- as.factor(matrix.long$value)

ggplot(matrix.long, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white") +
  theme_linedraw() +
  scale_fill_manual(values = c("white", "black"), name = "value") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  xlab("ASV") + ylab("Sample") +
  ggtitle("Consistently detected ASVs", 
          subtitle = "Determined through assessing only PCR reps that had >= 1% of total retained read depth (after running the pipeline)")

#histogram----
#plot
ggplot(all_1perc_seqs, aes(x=Sequence_length)) +
  geom_histogram(binwidth = 1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  xlab("Sequence length (in bp)") + ylab("Number of ASVs this length") +
  ggtitle("ASV sequence length distribution")



#stacked bar graph----
#plot total depth
ggplot(all_1perc_seqs, aes(x=Sample, y=Total.depth, fill=ASV)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("Sample") + ylab("Total read depth retained") +
  ggtitle("Proportion of ASVs within each sample")

#plot average depth
ggplot(all_1perc_seqs, aes(x=Sample, y=Average.depth, fill=ASV)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("Fish") + ylab("Average read depth retained") +
  ggtitle("Proportion of ASVs within each sample")

#get some data
total <- aggregate(Total.depth ~ ASV, data=all_1perc_seqs, FUN=sum)

#make fasta for downstream analysis----
#make data frame of only ASV and seq then convert to .fasta
#pull out only unique ASVs
true_seqs <- aggregate(Total.depth ~ ASV + Sequence, data = all_1perc_seqs, FUN = sum)

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
            file = "eDNA/aquatron/230509/ASVs.fasta", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)



#alig the sequences----
# load the sequences from the file
# change "DNA" to "RNA" or "AA" if necessary
library("Biostrings")
library("here")
#"here() starts at /Users/samanthabeal/Documents/MSc/Bioinformatics"

ASVs1 <- readDNAStringSet(here("eDNA/aquatron/230509/ASVs.fasta"), format = "fasta")

# look at some of the sequences (optional)
ASVs1

# nucleotide sequences need to be in the same orientation
library("DECIPHER")
ASVs1 <- OrientNucleotides(ASVs1)

# perform the alignment
alignedASVs <- AlignSeqs(ASVs1)

# view the alignment in a browser (optional)
BrowseSeqs(alignedASVs, highlight=0)

# write the alignment to a new FASTA file
writeXStringSet(alignedASVs,
                file="eDNA/aquatron/230509/alignedASVs.fasta")

# Find differences among sequences (Hamming distance)----
library(matrixcalc)
ASVs1_diffs <- as.matrix(stringDist(alignedASVs,method="hamming"))

# Visualize distance matrix

#Concert matrix to long form to work with ggplot heatmap
library(tibble)
ASVs1_diffs_long <- ASVs1_diffs %>% as.data.frame() %>%
  rownames_to_column("ASV") %>%
  pivot_longer(-c(ASV), names_to="ASV2",values_to="distance")

ggplot(ASVs1_diffs_long, aes(x=ASV,y=ASV2,fill=distance))+
  geom_raster()+
  scale_fill_gradient(low="white",high="darkred") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("ASV1") + ylab("ASV2") +
  ggtitle("Base differences among UNIX SmaI-corII ASVs")



#make haplotype network----
library("ape")
library("pegas")

UNIXdata <-read.dna("230508/alignedASVs.fasta", format="fasta")
UNIXhaplo <- haplotype(UNIXdata)
UNIXhaplo

# Number of haplotypes:5, sequence length-92

UNIXhaplodist <- dist.dna(UNIXhaplo, "N")
UNIXnetwork <- rmst(UNIXhaplodist, quiet = TRUE)
UNIXnetwork
#5 haplotypes
#9 links
#link lengths between 0 and 1 steps

plot(UNIXnetwork)
plot(UNIXnetwork, fast = TRUE)

#use popart to make haplotype network: need to convert to nexus file
#make sure in ~samanthabeal directory in terminal
#seqmagick convert --output-format nexus --alphabet dna Documents/MSc/Bioinformatics/UNIX-opt-new/eDNA/aquatron/230509/alignedASVs.fasta Documents/MSc/Bioinformatics/UNIX-opt-new/eDNA/aquatron/230509/alignedASVs.nex


