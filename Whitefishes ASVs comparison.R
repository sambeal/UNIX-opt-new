#compare ASVs from different sources
#gDNA = atlantic whitefish
#aquatron and release = atlantic whitefish
#shingle = lake whitefish

#Aquatron----
#read in data
aqua <- read.csv("eDNA/aquatron/230509/1%_ASVs.csv")
aqua$ASV <- as.factor(aqua$ASV)
aqua$Sequence <- as.factor(aqua$Sequence)
aqua$Sample <- as.factor(aqua$Sample)

#rename ASV so can add to the same 
aqua$ASV <- paste0(aqua$ASV, "-AW")

#plot
library(ggplot2)
#histogram
ggplot(aqua, aes(x=Sequence_length)) +
  geom_histogram(binwidth = 1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  xlab("Sequence length (in bp)") + ylab("Number of ASVs this length") +
  ggtitle("Aquatron ASV sequence length distribution")



#stacked bar graph
#plot total depth per sample
ggplot(aqua, aes(x=Sample, y=Total.depth, fill=ASV)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("Sample") + ylab("Total read depth retained") +
  ggtitle("Proportion of ASVs within each Aquatron sample")

#plot total depth across all samples
ggplot(aqua, aes(x="", y=Total.depth, fill=ASV)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("ASP02_02") + ylab("Total read depth retained") +
  ggtitle("Proportion of ASVs across all Aquatron samples")

#add percentage labels
aqua_percentage <- aggregate(Total.depth ~ ASV, data=aqua, FUN=sum)
aqua_percentage$Percent <- (aqua_percentage$Total.depth/sum(aqua_percentage$Total.depth))*100

ggplot(aqua_percentage, aes(x="", y=Total.depth, fill=ASV)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("ASP02_02") + ylab("Total read depth retained") +
  ggtitle("Proportion of ASVs across all Aquatron samples") +
  geom_text(aes(label = round(Percent, digits = 1)), size = 3, position = position_stack(vjust = 0.5))


#Shingle----
#read in data
shingle <- read.csv("eDNA/shingle/230510/shingleonly1%_ASVs.csv")
shingle$ASV <- as.factor(shingle$ASV)
shingle$Sequence <- as.factor(shingle$Sequence)
shingle$Sample <- as.factor(shingle$Sample)

#rename ASV so can add to the same 
shingle$ASV <- paste0(shingle$ASV, "-LW")

#plot
library(ggplot2)
#histogram
ggplot(shingle, aes(x=Sequence_length)) +
  geom_histogram(binwidth = 1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  xlab("Sequence length (in bp)") + ylab("Number of ASVs this length") +
  ggtitle("Shingle Lake ASV sequence length distribution")


#stacked bar graph
#plot total depth
#stacked bar graph
#plot total depth per sample
ggplot(shingle, aes(x=Sample, y=Total.depth, fill=ASV)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("Sample") + ylab("Total read depth retained") +
  ggtitle("Proportion of ASVs within each Shingle sample")

#plot total depth across all samples
ggplot(shingle, aes(x="", y=Total.depth, fill=ASV)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("Shingle") + ylab("Total read depth retained") +
  ggtitle("Proportion of ASVs across all Shingle samples")

#add percentage labels
shingle_percentage <- aggregate(Total.depth ~ ASV, data=shingle, FUN=sum)
shingle_percentage$Percent <- (shingle_percentage$Total.depth/sum(shingle_percentage$Total.depth))*100

ggplot(shingle_percentage, aes(x="", y=Total.depth, fill=ASV)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("Shingle") + ylab("Total read depth retained") +
  ggtitle("Proportion of ASVs across all Shingle samples") +
  geom_text(aes(label = round(Percent, digits = 1)), size = 3, position = position_stack(vjust = 0.5))


#Aquatron v. Shingle Comparison----
#want to compare them on a single graph
#need to make each unique sequence have its own name so they can be directly compared
#aqua
aqua_comp <- aqua_percentage
aqua_comp$Sample <- "ASP02_02"
aqua_seqs <- aqua[c(1:2)]
aqua_seqs <- unique(aqua_seqs)

library(tidyverse)
aqua_comp <- aqua_comp %>% 
  inner_join(aqua_seqs, by = c("ASV" = "ASV"))

#shingle
shingle_comp <- shingle_percentage
shingle_comp$Sample <- "Shingle"
shingle_seqs <- shingle[c(1:2)]
shingle_seqs <- unique(shingle_seqs)

library(tidyverse)
shingle_comp <- shingle_comp %>% 
  inner_join(shingle_seqs, by = c("ASV" = "ASV"))

#join
comp <- rbind(aqua_comp,shingle_comp)
comp$ASV <- as.factor(comp$ASV)
comp$Sequence <- as.factor(comp$Sequence)
comp$Sample <- as.factor(comp$Sample)

#subset into a dataframe for each sequence so can rename to be the same?
comp_split <- split(comp, comp$Sequence)
str(comp_split)
ASV_05 <- comp_split[[1]]
ASV_05$Name <- "ASV_05"

ASV_07 <- comp_split[[2]]
ASV_07$Name <- "ASV_07"

ASV_06 <- comp_split[[3]]
ASV_06$Name <- "ASV_06"

ASV_04 <- comp_split[[4]]
ASV_04$Name <- "ASV_04"

ASV_03 <- comp_split[[5]]
ASV_03$Name <- "ASV_03"

ASV_02 <- comp_split[[6]]
ASV_02$Name <- "ASV_02"

ASV_11 <- comp_split[[7]]
ASV_11$Name <- "ASV_11"

ASV_01 <- comp_split[[8]]
ASV_01$Name <- "ASV_01"

ASV_10 <- comp_split[[9]]
ASV_10$Name <- "ASV_10"

ASV_08 <- comp_split[[10]]
ASV_08$Name <- "ASV_08"

ASV_09 <- comp_split[[11]]
ASV_09$Name <- "ASV_09"

comp_renamed <- rbind(ASV_01, ASV_02, ASV_03, ASV_04,
                      ASV_05, ASV_06, ASV_07, ASV_08,
                      ASV_09, ASV_10, ASV_11)
comp_renamed$Name <- as.factor(comp_renamed$Name)

#save
library(readr)
write_csv(comp_renamed, "Aquatron and Shingle ASVs comparison.csv")

#plot stacked bar
ggplot(comp_renamed, aes(fill=Name, y=Total.depth, x=Sample)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("Sample") + ylab("Total read depth retained") +
  ggtitle("Proportion of ASVs across all each sample type") +
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

            

#Release----
#read in data
release <- read.csv("eDNA/release_02/230512/releaseonly1%_ASVs.csv")
release$ASV <- as.factor(release$ASV)
release$Sequence <- as.factor(release$Sequence)
release$Sample <- as.factor(release$Sample)

#rename ASV so can add to the same 
release$ASV <- paste0(release$ASV, "-AW")

#plot
library(ggplot2)
#histogram
ggplot(release, aes(x=Sequence_length)) +
  geom_histogram(binwidth = 1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  xlab("Sequence length (in bp)") + ylab("Number of ASVs this length") +
  ggtitle("Net pen release (wk2) ASV sequence length distribution")


#stacked bar graph
#plot total depth
#stacked bar graph
#plot total depth per sample
ggplot(release, aes(x=Sample, y=Total.depth, fill=ASV)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("Sample") + ylab("Total read depth retained") +
  ggtitle("Proportion of ASVs within each release sample")

#plot total depth across all samples
ggplot(release, aes(x="", y=Total.depth, fill=ASV)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("Shingle") + ylab("Total read depth retained") +
  ggtitle("Proportion of ASVs across all release samples")

#add percentage labels
release_percentage <- aggregate(Total.depth ~ ASV, data=release, FUN=sum)
release_percentage$Percent <- (release_percentage$Total.depth/sum(release_percentage$Total.depth))*100

ggplot(release_percentage, aes(x="", y=Total.depth, fill=ASV)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("Shingle") + ylab("Total read depth retained") +
  ggtitle("Proportion of ASVs across all release samples") +
  geom_text(aes(label = round(Percent, digits = 1)), size = 3, position = position_stack(vjust = 0.5))


#Aquatron, Shingle, and Release Comparison----
#want to compare them on a single graph
#need to make each unique sequence have its own name so they can be directly compared
#aqua
aqua_comp <- aqua_percentage
aqua_comp$Sample <- "Aquatron"
aqua_seqs <- aqua[c(1:2)]
aqua_seqs <- unique(aqua_seqs)

#library(tidyverse)
aqua_comp <- aqua_comp %>% 
  inner_join(aqua_seqs, by = c("ASV" = "ASV"))

#shingle
shingle_comp <- shingle_percentage
shingle_comp$Sample <- "Shingle"
shingle_seqs <- shingle[c(1:2)]
shingle_seqs <- unique(shingle_seqs)

#library(tidyverse)
shingle_comp <- shingle_comp %>% 
  inner_join(shingle_seqs, by = c("ASV" = "ASV"))

#release
release_comp <- release_percentage
release_comp$Sample <- "Release"
release_seqs <- release[c(1:2)]
release_seqs <- unique(release_seqs)

library(tidyverse)
release_comp <- release_comp %>% 
  inner_join(release_seqs, by = c("ASV" = "ASV"))

#join
comp <- rbind(aqua_comp,shingle_comp,release_comp)
comp$ASV <- as.factor(comp$ASV)
comp$Sequence <- as.factor(comp$Sequence)
comp$Sample <- as.factor(comp$Sample)

#subset into a dataframe for each sequence so can rename to be the same
comp_split <- split(comp, comp$Sequence)
str(comp_split)
ASV_05 <- comp_split[[1]]
ASV_05$Name <- "ASV_05"

ASV_06 <- comp_split[[2]]
ASV_06$Name <- "ASV_06"

ASV_04 <- comp_split[[3]]
ASV_04$Name <- "ASV_04"

ASV_03 <- comp_split[[4]]
ASV_03$Name <- "ASV_03"

ASV_02 <- comp_split[[5]]
ASV_02$Name <- "ASV_02"

ASV_09 <- comp_split[[6]]
ASV_09$Name <- "ASV_09"

ASV_01 <- comp_split[[7]]
ASV_01$Name <- "ASV_01"

ASV_08 <- comp_split[[8]]
ASV_08$Name <- "ASV_08"

ASV_07 <- comp_split[[9]]
ASV_07$Name <- "ASV_07"

ASV_10 <- comp_split[[10]]
ASV_10$Name <- "ASV_10"

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
write_csv(comp_renamed, "Aquatron, Shingle, and Release ASVs comparison.csv")

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
ggplot(comp_renamed, aes(x=Name, y=Total.depth)) +
  geom_bar(stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("Sample") + ylab("Total read depth retained") +
  ggtitle("Read depth of ASVs across all each sample type") +
  facet_grid(Sample~.)


aggregate(Total.depth ~ Name + Sample, data = comp_renamed, FUN = sum)

#faceted bar by percentage
aqua_percentage$Sample <- "Aquatron"
aqua_percentage$Species <- "Atlantic Whitefish"

shingle_percentage$Sample <- "Shingle"
shingle_percentage$Species <- "Lake Whitefish"

release_percentage$Sample <- "Release"
release_percentage$Species <- "Atlantic Whitefish"

comp_percent <- rbind(aqua_percentage, shingle_percentage, release_percentage)

comp_ASVs <- comp_renamed[c(1,6)]

comp_percent <- comp_percent %>% 
  inner_join(comp_ASVs, by = c("ASV" = "ASV"))

ggplot(comp_percent, aes(x=Name, y=Percent, fill=Species)) +
  geom_bar(stat = "identity") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) +
  xlab("Sample") + ylab("Percentage of read depth") +
  ggtitle("Read depth of ASVs across all each sample type") +
  facet_grid(Sample~.)


