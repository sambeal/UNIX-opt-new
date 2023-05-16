#analysis of 230508 ASVs 
all_1perc_seqs <- read.csv("230508/1%_ASVs.csv")

#plotting----
#fish/ASV tile matrix----
#pull out only fish and ASV
fish_ASVs <- all_1perc_seqs[c(1,7)]
fish_ASVs$ASV <- as.factor(fish_ASVs$ASV)
fish_ASVs$Fish <- as.factor(fish_ASVs$Fish)

uniqueASV <- as.data.frame(unique(fish_ASVs$ASV))
names(uniqueASV)[1] <- "ASV"
#6 unique ASVs

#make dataframe of unique fish for the matrix plot
uniqueFish <- as.data.frame(unique(fish_ASVs$Fish))
names(uniqueFish)[1] <- "Fish"

#subset into a dataframe for each ASV (less ASVs than fish)
ASV_split <- split(all_1perc_seqs, all_1perc_seqs$ASV)
str(ASV_split)
uniq1 <- ASV_split[[1]]
uniq2 <- ASV_split[[2]]
uniq3 <- ASV_split[[3]]
uniq4 <- ASV_split[[4]]
uniq5 <- ASV_split[[5]]
uniq6 <- ASV_split[[6]]


#check if fish in in each ASV df
#uniq1
uniq1_fish <- as.data.frame(uniqueFish$Fish %in% uniq1$Fish)
names(uniq1_fish)[1] <- "uniq1"

#uniq2
uniq2_fish <- as.data.frame(uniqueFish$Fish %in% uniq2$Fish)
names(uniq2_fish)[1] <- "uniq2"

#uniq3
uniq3_fish <- as.data.frame(uniqueFish$Fish %in% uniq3$Fish)
names(uniq3_fish)[1] <- "uniq3"

#uniq4
uniq4_fish <- as.data.frame(uniqueFish$Fish %in% uniq4$Fish)
names(uniq4_fish)[1] <- "uniq4"

#uniq5
uniq5_fish <- as.data.frame(uniqueFish$Fish %in% uniq5$Fish)
names(uniq5_fish)[1] <- "uniq5"

#uniq6
uniq6_fish <- as.data.frame(uniqueFish$Fish %in% uniq6$Fish)
names(uniq6_fish)[1] <- "uniq6"

#add together
ASVperfish <- cbind(uniq1_fish,uniq2_fish,uniq3_fish,
                    uniq4_fish, uniq5_fish, uniq6_fish)
row.names(ASVperfish) <- uniqueFish$Fish

#save
library(readr)
write_csv(ASVperfish, "230508/true_false_matrix.csv")

#change from true/false to 1/0 by making a matrix where true=1, false=0
ASVperfishmatrix <- data.matrix(ASVperfish)

#plot
library(ggplot2)
library(reshape2)

#Convert wide to long
matrix.long <- melt(ASVperfishmatrix)
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
  xlab("ASV") + ylab("Fish") +
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
ggplot(all_1perc_seqs, aes(x=Fish, y=Total.depth, fill=ASV)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle=70, vjust = 0.6)) +
  xlab("Fish") + ylab("Total read depth retained") +
  ggtitle("Proportion of ASVs within each fish")

#plot average depth
ggplot(all_1perc_seqs, aes(x=Fish, y=Average.depth, fill=ASV)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle=70, vjust = 0.6)) +
  xlab("Fish") + ylab("Average read depth retained") +
  ggtitle("Proportion of ASVs within each fish")



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
            file = "230508/ASVs.fasta", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)



#alig the sequences----
# load the sequences from the file
# change "DNA" to "RNA" or "AA" if necessary
library("Biostrings")
library("here")
#"here() starts at /Users/samanthabeal/Documents/MSc/Bioinformatics"

ASVs1 <- readDNAStringSet(here("230508/ASVs.fasta"), format = "fasta")

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
                file="230508/alignedASVs.fasta")

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
#seqmagick convert --output-format nexus --alphabet dna Documents/MSc/Bioinformatics/UNIX-opt-new/230508/alignedASVs.fasta Documents/MSc/Bioinformatics/UNIX-opt-new/230508/alignedASVs.nex


