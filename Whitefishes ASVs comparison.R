#compare ASVs from different sources
#gDNA = atlantic whitefish
#aquatron and release = atlantic whitefish
#shingle = lake whitefish

#working directory = "/Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-opt-new"

#gDNA----
#read in data
gDNA <- read.csv("gDNA/230516/1%_ASVs.csv")
gDNA$Sample <- "Atlantic Whitefish gDNA"

gDNA$ASV <- as.factor(gDNA$ASV)
gDNA$Sequence <- as.factor(gDNA$Sequence)
gDNA$Sample <- as.factor(gDNA$Sample)

#rename ASV so can add to the same df with other samples later
gDNA$ASV <- paste0(gDNA$ASV, "-AW")

#plot
library(ggplot2)
gDNAuniq <- gDNA[c(1:3)]
gDNAuniq <- as.data.frame(unique(gDNAuniq))

#histogram
ggplot(gDNAuniq, aes(x=Sequence_length)) +
  geom_histogram(binwidth = 1, colour = "grey", fill = "black") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  xlab("Sequence length (bp)") + ylab("Number of ASVs this length") +
  ggtitle("gDNA ASV sequence length distribution")

#stacked bar graph
#plot total depth as % of reads per fish
ggplot(gDNA, aes(x=Fish, y=Total.depth, fill=ASV)) +
  geom_bar(position = "fill", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle=70, vjust = 0.6)) +
  xlab("Fish") + ylab("Proportion of retained depth") +
  ggtitle("Proportion of ASVs within each fish")



ggplot(gDNA, aes(x="", y=Total.depth, fill=ASV)) +
  geom_bar(position = "fill", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("gDNA from 16 Atlantic Whitefish") + ylab("Proportion of retained depth") +
  ggtitle("Proportion of ASVs across all gDNA samples") 

#Aquatron----
#read in data
aqua <- read.csv("eDNA/aquatron/230516/1%_ASVs.csv")
names(aqua)[7] <- "Sample_Name"
aqua$Sample <- "Aquatron"

aqua$ASV <- as.factor(aqua$ASV)
aqua$Sequence <- as.factor(aqua$Sequence)
aqua$Sample_Name <- as.factor(aqua$Sample_Name)
aqua$Sample <- as.factor(aqua$Sample)


#rename ASV so can add to the same 
aqua$ASV <- paste0(aqua$ASV, "-AW") 

#plot
library(ggplot2)
#histogram
aquauniq <- aqua[c(1:3)]
aquauniq <- as.data.frame(unique(aquauniq))

#histogram
ggplot(aquauniq, aes(x=Sequence_length)) +
  geom_histogram(binwidth = 1, colour = "grey", fill = "black") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  xlab("Sequence length (bp)") + ylab("Number of ASVs this length") +
  ggtitle("Aquatron ASV sequence length distribution")


#stacked bar graph
#plot total depth per sample
ggplot(aqua, aes(x=Sample_Name, y=Total.depth, fill=ASV)) +
  geom_bar(position = "fill", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("Sample") + ylab("Proportion of retained depth") +
  ggtitle("Proportion of ASVs within each Aquatron sample")

#plot total depth across all samples
ggplot(aqua, aes(x="", y=Total.depth, fill=ASV)) +
  geom_bar(position = "fill", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("ASP02_02") + ylab("Proportion of retained depth") +
  ggtitle("Proportion of ASVs across all Aquatron samples")

#to add percentage labels = need to make new df, use above code to plot
#but change data from aqua to aqua_percentage when plotting
#aqua_percentage <- aggregate(Total.depth ~ ASV, data=aqua, FUN=sum)
#aqua_percentage$Percent <- (aqua_percentage$Total.depth/sum(aqua_percentage$Total.depth))*100
#+ geom_text(aes(label = round(Percent, digits = 1)), size = 3, position = position_stack(vjust = 0.5))


#Shingle----
#read in data
shingle <- read.csv("eDNA/shingle/230516/1%_ASVs.csv")
names(shingle)[7] <- "Sample_Name"
shingle$Sample <- "Shingle Lake"
shingle$ASV <- as.factor(shingle$ASV)
shingle$Sequence <- as.factor(shingle$Sequence)
shingle$Sample <- as.factor(shingle$Sample)
shingle$Sample_Name <- as.factor(shingle$Sample_Name)

#rename ASV so can add to the same 
shingle$ASV <- paste0(shingle$ASV, "-LW")

#plot
library(ggplot2)
#histogram
shingleuniq <- shingle[c(1:3)]
shingleuniq <- as.data.frame(unique(shingleuniq))

#histogram
ggplot(shingleuniq, aes(x=Sequence_length)) +
  geom_histogram(binwidth = 1, colour = "grey", fill = "black") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  xlab("Sequence length (bp)") + ylab("Number of ASVs this length") +
  ggtitle("Shingle ASV sequence length distribution")


#stacked bar graph
#plot total depth per sample
ggplot(shingle, aes(x=Sample_Name, y=Total.depth, fill=ASV)) +
  geom_bar(position = "fill", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("Sample") + ylab("Proportion of retained depth") +
  ggtitle("Proportion of ASVs within each Shingle sample")

#plot total depth across all samples
ggplot(shingle, aes(x="", y=Total.depth, fill=ASV)) +
  geom_bar(position = "fill", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("Shingle Lake") + ylab("Proportion of retained depth") +
  ggtitle("Proportion of ASVs across all Shingle samples")


#Release----
#read in data
release <- read.csv("eDNA/release_02/230516/1%_ASVs.csv")
names(release)[7] <- "Sample_Name"
release$Sample <- "Net pen release"
release$ASV <- as.factor(release$ASV)
release$Sequence <- as.factor(release$Sequence)
release$Sample <- as.factor(release$Sample)
release$Sample_Name <- as.factor(release$Sample_Name)

#rename ASV so can add to the same 
release$ASV <- paste0(release$ASV, "-AW")

#plot
library(ggplot2)
#histogram
releaseuniq <- release[c(1:3)]
releaseuniq <- as.data.frame(unique(releaseuniq))

#histogram
ggplot(releaseuniq, aes(x=Sequence_length)) +
  geom_histogram(binwidth = 1, colour = "grey", fill = "black") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  xlab("Sequence length (bp)") + ylab("Number of ASVs this length") +
  ggtitle("Release ASV sequence length distribution")


#stacked bar graph
#plot total depth per sample
ggplot(release, aes(x=Sample, y=Total.depth, fill=ASV)) +
  geom_bar(position = "fill", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("Sample") + ylab("Proportion of retained depth") +
  ggtitle("Proportion of ASVs within each Release sample")

#plot total depth across all samples
ggplot(release, aes(x="", y=Total.depth, fill=ASV)) +
  geom_bar(position = "fill", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("Release") + ylab("Proportion of retained depth") +
  ggtitle("Proportion of ASVs across all Release samples")

#gDNA, Aquatron, Shingle, and Release Comparison----
#want to compare them on a single graph
#need to make each unique sequence have its own name so they can be directly compared
#gDNA
gDNA_percentage <- aggregate(Total.depth ~ ASV, data=gDNA, FUN=sum)
gDNA_percentage$Percent <- (gDNA_percentage$Total.depth/sum(gDNA_percentage$Total.depth))*100

gDNA_comp <- gDNA_percentage
gDNA_comp$Sample <- "Atlantic Whitefish gDNA"
gDNA_seqs <- gDNA[c(1:2)]
gDNA_seqs <- unique(gDNA_seqs)

library(tidyverse)
gDNA_comp <- gDNA_comp %>% 
  inner_join(gDNA_seqs, by = c("ASV" = "ASV"))

#aqua
aqua_percentage <- aggregate(Total.depth ~ ASV, data=aqua, FUN=sum)
aqua_percentage$Percent <- (aqua_percentage$Total.depth/sum(aqua_percentage$Total.depth))*100

aqua_comp <- aqua_percentage
aqua_comp$Sample <- "Aquatron"
aqua_seqs <- aqua[c(1:2)]
aqua_seqs <- unique(aqua_seqs)

library(tidyverse)
aqua_comp <- aqua_comp %>% 
  inner_join(aqua_seqs, by = c("ASV" = "ASV"))

#shingle
shingle_percentage <- aggregate(Total.depth ~ ASV, data=shingle, FUN=sum)
shingle_percentage$Percent <- (shingle_percentage$Total.depth/sum(shingle_percentage$Total.depth))*100

shingle_comp <- shingle_percentage
shingle_comp$Sample <- "Shingle"
shingle_seqs <- shingle[c(1:2)]
shingle_seqs <- unique(shingle_seqs)

library(tidyverse)
shingle_comp <- shingle_comp %>% 
  inner_join(shingle_seqs, by = c("ASV" = "ASV"))

#release
release_percentage <- aggregate(Total.depth ~ ASV, data=release, FUN=sum)
release_percentage$Percent <- (release_percentage$Total.depth/sum(release_percentage$Total.depth))*100

release_comp <- release_percentage
release_comp$Sample <- "Net pen release"
release_seqs <- release[c(1:2)]
release_seqs <- unique(release_seqs)

library(tidyverse)
release_comp <- release_comp %>% 
  inner_join(release_seqs, by = c("ASV" = "ASV"))

#join
comp <- rbind(gDNA_comp, aqua_comp,shingle_comp,release_comp)
comp$ASV <- as.factor(comp$ASV)
comp$Sequence <- as.factor(comp$Sequence)
comp$Sample <- as.factor(comp$Sample)

#subset into a dataframe for each sequence so can rename to be the same
comp_split <- split(comp, comp$Sequence)
str(comp_split)
ASV_05 <- comp_split[[1]]
ASV_05$Name <- "ASV_05"

ASV_04 <- comp_split[[2]]
ASV_04$Name <- "ASV_04"

ASV_03 <- comp_split[[3]]
ASV_03$Name <- "ASV_03"

ASV_02 <- comp_split[[4]]
ASV_02$Name <- "ASV_02"

ASV_10 <- comp_split[[5]]
ASV_10$Name <- "ASV_10"

ASV_01 <- comp_split[[6]]
ASV_01$Name <- "ASV_01"

ASV_06 <- comp_split[[7]]
ASV_06$Name <- "ASV_06"

ASV_08 <- comp_split[[8]]
ASV_08$Name <- "ASV_08"

ASV_09 <- comp_split[[9]]
ASV_09$Name <- "ASV_09"

ASV_07 <- comp_split[[10]]
ASV_07$Name <- "ASV_07"

ASV_11 <- comp_split[[11]]
ASV_11$Name <- "ASV_11"
  
ASV_12 <- comp_split[[12]]
ASV_12$Name <- "ASV_12"

comp_renamed <- rbind(ASV_01, ASV_02, ASV_03, ASV_04,
                      ASV_05, ASV_06, ASV_07, ASV_08,
                      ASV_09, ASV_10, ASV_11, ASV_12)

#add back in seq lengths
comp_renamed$Sequence_length <- str_length(comp_renamed$Sequence)

comp_renamed$Name <- as.factor(comp_renamed$Name)

#save
library(readr)
write_csv(comp_renamed, "gDNA, Aquatron, Shingle, and Release ASVs comparison.csv")

#plot stacked bar
ggplot(comp_renamed, aes(fill=Name, y=Total.depth, x=Sample)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("Sample") + ylab("Total read depth retained") +
  ggtitle("Read depth of ASVs across all each sample type") +
  geom_text(aes(label = round(Percent, digits = 1)), size = 3, position = position_stack(vjust = 0.5))

#percentage bar
ggplot(comp_renamed, aes(fill=Name, y=Total.depth, x=Sample)) +
  geom_bar(position = "fill", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("Sample") + ylab("Percentage of total read depth") +
  ggtitle("Proportion of ASVs across all each sample type")

#faceted bar
ggplot(comp_renamed, aes(x=Name, y=Percent)) +
  geom_bar(stat = "identity", fill = "black") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("ASV") + ylab("Percent of depth retained") +
  ggtitle("ASV proportions across each sample type") +
  facet_grid(Sample~.)

#add species and colour by that
gDNA_percentage$Sample <- "Atlantic Whitefish gDNA"
gDNA_percentage$Species <- "Atlantic Whitefish"

aqua_percentage$Sample <- "Aquatron"
aqua_percentage$Species <- "Atlantic Whitefish"

shingle_percentage$Sample <- "Shingle Lake"
shingle_percentage$Species <- "Lake Whitefish"

release_percentage$Sample <- "Atlantic Whitefish Release"
release_percentage$Species <- "Atlantic Whitefish"

comp_percent <- rbind(gDNA_percentage, aqua_percentage, shingle_percentage, release_percentage)

comp_ASVs <- comp_renamed[c(1,6)]

comp_percent <- comp_percent %>% 
  inner_join(comp_ASVs, by = c("ASV" = "ASV"))

ggplot(comp_percent, aes(x=Name, y=Percent, fill = Species)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("ASV") + ylab("Percent of depth retained") +
  ggtitle("ASV proportions across each sample type") +
  facet_grid(Sample~.)

#make a fasta of just these sequences for later analysis
fasta <- unique(comp_renamed[5:6])
library(tidyverse)
fasta <- fasta %>% relocate(Sequence, .after = Name)

fasta$Name <- paste0(">", fasta$Name)

library(tidyverse)

fasta_print <-
  fasta %>%
  select(Name,Sequence) %>%
  rowwise() %>%
  pivot_longer(Name:Sequence) %>%
  select(-name)

#allseqs
write.table(fasta_print, 
            file = "whitefish_comparison_ASVs.fasta", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)


#230613----
#need to rename the gDNA ASVs to match all others (aka uniq1 = ASV_01)
#compile all AW ASVs into one facet, LW into another (for better simplicity)

renamedASVs <- comp_renamed[c(5:6)]
renamedASVsunique <- as.data.frame(unique(renamedASVs))

#gDNA
gDNAnew <- gDNA
gDNAnew <- gDNAnew %>% 
  inner_join(renamedASVsunique, by = c("Sequence" = "Sequence"))
names(gDNAnew)[7] <- "Sample_Name"

gDNAnew$Name <- gsub("ASV_", "", gDNAnew$Name)
gDNAnew$Sample_Name <- gsub("1909", "", gDNAnew$Sample_Name)
gDNAnew$Sample_Name <- gsub("1902", "", gDNAnew$Sample_Name)



#stacked bar graph
#plot total depth as % of reads per fish
ggplot(gDNAnew, aes(x=Sample_Name, y=Total.depth, fill=Name)) +
  geom_bar(position = "fill", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20),
        axis.text.x = element_text(vjust = 0.6)) +
  xlab("Fish") + ylab("Proportion of retained depth") +
  labs(fill = "ASV")

#aqua
aquanew <- aqua
aquanew <- aquanew %>% 
  inner_join(renamedASVsunique, by = c("Sequence" = "Sequence"))

#stacked bar graph
#plot total depth as % of reads per fish
ggplot(aquanew, aes(x=Sample_Name, y=Total.depth, fill=Name)) +
  geom_bar(position = "fill", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("Sample") + ylab("Proportion of retained depth") +
  labs(fill = "ASV")

#shingle
shinglenew <- shingle
shinglenew <- shinglenew %>% 
  inner_join(renamedASVsunique, by = c("Sequence" = "Sequence"))

#release
releasenew <- release
releasenew <- releasenew %>% 
  inner_join(renamedASVsunique, by = c("Sequence" = "Sequence"))

ggplot(releasenew, aes(x=Sample_Name, y=Total.depth, fill=Name)) +
  geom_bar(position = "fill", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("Sample") + ylab("Proportion of retained depth") +
  labs(fill = "ASV")


#join AWF eDNA samples only
awfeDNA <- rbind(aquanew, releasenew)
awfeDNA$Name <- gsub("ASV_", "", awfeDNA$Name)



ggplot(awfeDNA, aes(x=Sample_Name, y=Total.depth, fill=Name)) +
  geom_bar(position = "fill", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("Sample") + ylab("Proportion of retained depth") +
  labs(fill = "ASV")

#make fasta so can get alignment
awfeDNAseqs <- unique(awfeDNA[c(2,12)])
library(dplyr)
awfeDNAseqs <- awfeDNAseqs %>% relocate(Name, .before = Sequence)
awfeDNAseqs$Name <- paste0(">", awfeDNAseqs$Name)

awfeDNAfasta <-
  awfeDNAseqs %>%
  select(Name,Sequence) %>%
  rowwise() %>%
  pivot_longer(Name:Sequence) %>%
  select(-name)

#allseqs
write.table(awfeDNAfasta, 
            file = "eDNA/230614_AWF_eDNA_ASVs.fasta", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)


#add species name and join
#add species and colour by that
gDNAnew$Species <- "Atlantic Whitefish"

aquanew$Species <- "Atlantic Whitefish"

shinglenew$Species <- "Lake Whitefish"

releasenew$Species <- "Atlantic Whitefish"

newcomp <- rbind(gDNAnew, aquanew,shinglenew,releasenew)
newcomp$ASV <- as.factor(newcomp$ASV)
newcomp$Sequence <- as.factor(newcomp$Sequence)
newcomp$Sample <- as.factor(newcomp$Sample)
newcomp$Species <- as.factor(newcomp$Species)

#remove 'ASV_" from name to make it shorter and nicer to plot
newcomp$Name <- gsub("ASV_", "", newcomp$Name)



ggplot(newcomp, aes(x=Name, y=Total.depth, fill = Sample)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 30),
        axis.text.x = element_text(vjust = 0.6),
        strip.text = element_text(size = 30)) +
  xlab("ASV") + ylab("Total depth retained") +
  facet_grid(Species~.)

