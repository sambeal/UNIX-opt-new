#analysis of the ASVs returned from the UNIX pipeline when run with default parameters except:
#sequence length restricted to 70-105bp (as determined from multiqc report)
#phred score 26 enforced during trimming
#sequences only taken through the dereplication stage

#when working with field samples that included an aquatron positive control, remove it prior to analysis
#(depths are higher for pos con and therefore will skew the 1% cut off - want to keep samples that are 
# >=1% of the field reads)


#230508 1% ASV determination of gDNA
#assess holistically
#only keep samples (individual PCR reps) that >=1% total retained read depth
#to remove low quality PCRs while maintaining high quality ones

#230509 1% ASV determination of eDNA (ASP02_02_001-3)

#230510 1% ASV determination of eDNA (SHI01_01_002-_02_003)

#230510 1% ASV determination of eDNA (AWR000_02_001-3)


#set working directory----
getwd()

#load data----
library(readr)
all <- readr::read_tsv("eDNA/release_02/output/ASV/derep.tsv", show_col_types = FALSE)
all <- as.data.frame(all)

#rename cols to be only fish name + index combo
names(all) = gsub(pattern = "_.*", replacement = "", x = names(all))
names(all)[1] <- "ASV"

#it's being weird when the ASV column is still there so will remove it and add it back before subsetting/fish
ASV <- all[c(1)]
all <- all[-c(1)]

#shingle+aqautron pos con----
#find total depth across all samples
total <- sum(all[1:19])

#cutoff=
total*0.01

#keep only samples that are >=1% of total depth
all_1perc <- all[colSums(all[1:19]) >= (total*0.01)]

#add back ASV column
all_1perc <- cbind(ASV, all_1perc)


#pull out individual samples
#ASP02_02_002 (pos con)
#2/3 PCR reps only
ASP02_02_002 <- all_1perc[c(1:3)]
names(ASP02_02_002)[2] <- "PCR.1"
names(ASP02_02_002)[3] <- "PCR.2"

#add column of sample ID
ASP02_02_002$Sample <- "ASP02_02_002"

#make new df of the 1% ASVs in this sample PER PCR REP (not across all reps)
ASP02_02_002_1perc <- subset(ASP02_02_002, PCR.1>=(sum(ASP02_02_002$PCR.1)*0.01) | PCR.2>=(sum(ASP02_02_002$PCR.2)*0.01))

#add empty PCR.3 so can bind all dfs together later
ASP02_02_002_1perc$PCR.3 <- NA

#reorganize to match other fish dfs
library(dplyr)
ASP02_02_002_1perc <- ASP02_02_002_1perc %>% relocate(PCR.3, .after = PCR.2)

#add total depth and average depth columns
ASP02_02_002_1perc$Total.depth <-apply(ASP02_02_002_1perc[,2:3],1,sum)
ASP02_02_002_1perc$Average.depth <-apply(ASP02_02_002_1perc[,2:3],1,mean)

\
#SHI01_01_002
SHI01_01_002 <- all_1perc[c(1,4:5)]
names(SHI01_01_002)[2] <- "PCR.1"
names(SHI01_01_002)[3] <- "PCR.2"

#add column of fish ID
SHI01_01_002$Sample <- "SHI01_01_002"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
SHI01_01_002_1perc <- subset(SHI01_01_002, PCR.1>=(sum(SHI01_01_002$PCR.1)*0.01) | PCR.2>=(sum(SHI01_01_002$PCR.2)*0.01))

#add empty PCR.3 so can bind all dfs together later
SHI01_01_002_1perc$PCR.3 <- NA

#reorganize to match other fish dfs
SHI01_01_002_1perc <- SHI01_01_002_1perc %>% relocate(PCR.3, .after = PCR.2)

#add total depth and average depth columns
SHI01_01_002_1perc$Total.depth <-apply(SHI01_01_002_1perc[,2:3],1,sum)
SHI01_01_002_1perc$Average.depth <-apply(SHI01_01_002_1perc[,2:3],1,mean)

#SHI01_01_003
SHI01_01_003 <- all_1perc[c(1,6:7)]
names(SHI01_01_003)[2] <- "PCR.1"
names(SHI01_01_003)[3] <- "PCR.2"

#add column of fish ID
SHI01_01_003$Sample <- "SHI01_01_003"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
SHI01_01_003_1perc <- subset(SHI01_01_003, PCR.1>=(sum(SHI01_01_003$PCR.1)*0.01) | PCR.2>=(sum(SHI01_01_003$PCR.2)*0.01))

#add empty PCR.3 so can bind all dfs together later
SHI01_01_003_1perc$PCR.3 <- NA

#reorganize to match other fish dfs
SHI01_01_003_1perc <- SHI01_01_003_1perc %>% relocate(PCR.3, .after = PCR.2)

#add total depth and average depth columns
SHI01_01_003_1perc$Total.depth <-apply(SHI01_01_003_1perc[,2:3],1,sum)
SHI01_01_003_1perc$Average.depth <-apply(SHI01_01_003_1perc[,2:3],1,mean)


#SHI02_01_001
SHI02_01_001 <- all_1perc[c(1,8:10)]
names(SHI02_01_001)[2] <- "PCR.1"
names(SHI02_01_001)[3] <- "PCR.2"
names(SHI02_01_001)[4] <- "PCR.3"

#add column of fish ID
SHI02_01_001$Sample <- "SHI02_01_001"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
SHI02_01_001_1perc <- subset(SHI02_01_001, PCR.1>=(sum(SHI02_01_001$PCR.1)*0.01) | PCR.2>=(sum(SHI02_01_001$PCR.2)*0.01) | PCR.3>=(sum(SHI02_01_001$PCR.3)*0.01))

#add total depth and average depth columns
SHI02_01_001_1perc$Total.depth <-apply(SHI02_01_001_1perc[,2:4],1,sum)
SHI02_01_001_1perc$Average.depth <-apply(SHI02_01_001_1perc[,2:4],1,mean)

#SHI02_01_002
SHI02_01_002 <- all_1perc[c(1,11:13)]
names(SHI02_01_002)[2] <- "PCR.1"
names(SHI02_01_002)[3] <- "PCR.2"
names(SHI02_01_002)[4] <- "PCR.3"

#add column of fish ID
SHI02_01_002$Sample <- "SHI02_01_002"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
SHI02_01_002_1perc <- subset(SHI02_01_002, PCR.1>=(sum(SHI02_01_002$PCR.1)*0.01) | PCR.2>=(sum(SHI02_01_002$PCR.2)*0.01) | PCR.3>=(sum(SHI02_01_002$PCR.3)*0.01))

#add total depth and average depth columns
SHI02_01_002_1perc$Total.depth <-apply(SHI02_01_002_1perc[,2:4],1,sum)
SHI02_01_002_1perc$Average.depth <-apply(SHI02_01_002_1perc[,2:4],1,mean)

#SHI02_01_003
SHI02_01_003 <- all_1perc[c(1,14:16)]
names(SHI02_01_003)[2] <- "PCR.1"
names(SHI02_01_003)[3] <- "PCR.2"
names(SHI02_01_003)[4] <- "PCR.3"

#add column of fish ID
SHI02_01_003$Sample <- "SHI02_01_003"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
SHI02_01_003_1perc <- subset(SHI02_01_003, PCR.1>=(sum(SHI02_01_003$PCR.1)*0.01) | PCR.2>=(sum(SHI02_01_003$PCR.2)*0.01) | PCR.3>=(sum(SHI02_01_003$PCR.3)*0.01))

#add total depth and average depth columns
SHI02_01_003_1perc$Total.depth <-apply(SHI02_01_003_1perc[,2:4],1,sum)
SHI02_01_003_1perc$Average.depth <-apply(SHI02_01_003_1perc[,2:4],1,mean)

#put it all together
all_1perc <- rbind(ASP02_02_002_1perc, SHI01_01_002_1perc, SHI01_01_003_1perc, 
                   SHI02_01_001_1perc, SHI02_01_002_1perc, SHI02_01_003_1perc)




#add sequences - make it new from shingle folder
library (devtools)
library(tidyverse)
source_url("https://raw.githubusercontent.com/lrjoshi/FastaTabular/master/fasta_and_tabular.R")

FastaToTabular("eDNA/shingle/output/ASV/derep.fasta")
#The output will be stored as dna_table.csv in the current directory ("Bioinformatics")

#read in data
UNIXseqs <- read.csv("dna_table.csv")


#reformat, remove first column, name new 1st column to ASV and new 2nd to Sequence
#need the columns to match so can join to all_1perc 
UNIXseqs <- UNIXseqs[-1]
names(UNIXseqs)[1] <- "ASV"
names(UNIXseqs)[2] <- "Sequence"

#rename the contents of the ASV column to match the ASV column of all_1perc
#";.*" says: replace semicolon (;) and every character after that (.*), with nothing "".
UNIXseqs$ASV <- gsub(";.*", "", UNIXseqs$ASV)
#">" says: replace arrow (>) with nothing "".
UNIXseqs$ASV <- gsub(">", "", UNIXseqs$ASV)

#add data frames together
#library(tidyverse)
all_1perc_seqs <- all_1perc %>% 
  inner_join(UNIXseqs, by = c("ASV" = "ASV"))

#move sequence to the front
all_1perc_seqs <- all_1perc_seqs %>% relocate(Sequence, .after = ASV)

#add column of # PCR reps the reads came from
all_1perc_seqs$PCR.reps <- ifelse(all_1perc_seqs$PCR.1>0 & !is.na(all_1perc_seqs$PCR.3), 3, 2)


#add seq length
all_1perc_seqs$Sequence_length <- str_length(all_1perc_seqs$Sequence)

#reformat
all_1perc_seqs <- all_1perc_seqs %>% relocate(Sequence_length, .after = Sequence)
depths <- aggregate(Total.depth ~ ASV, data=all_1perc_seqs, FUN=sum)


#save output for later analysis
write_csv(all_1perc_seqs, "eDNA/shingle/230510/1%_ASVs.csv")



#shingle samples only----
#try again and remove ASP samples from the data frame (pos con to know if PCR worked)
shingle <- all[-c(1:3)]

total_shingle <- sum(shingle[1:16])

#keep only samples that are >=1% of total shingle depth
shingle_1perc <- shingle[colSums(shingle[1:16]) >= (total_shingle*0.01)]

#add back ASV column
shingle_1perc <- cbind(ASV, shingle_1perc)


#pull out individual samples
#SHI01_01_002
SHI01_01_002 <- shingle_1perc[c(1:3)]
names(SHI01_01_002)[2] <- "PCR.1"
names(SHI01_01_002)[3] <- "PCR.2"

#add column of fish ID
SHI01_01_002$Sample <- "SHI01_01_002"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
SHI01_01_002_1perc <- subset(SHI01_01_002, PCR.1>=(sum(SHI01_01_002$PCR.1)*0.01) | PCR.2>=(sum(SHI01_01_002$PCR.2)*0.01))

#add empty PCR.3 so can bind all dfs together later
SHI01_01_002_1perc$PCR.3 <- NA

#reorganize to match other fish dfs
SHI01_01_002_1perc <- SHI01_01_002_1perc %>% relocate(PCR.3, .after = PCR.2)

#add total depth and average depth columns
SHI01_01_002_1perc$Total.depth <-apply(SHI01_01_002_1perc[,2:3],1,sum)
SHI01_01_002_1perc$Average.depth <-apply(SHI01_01_002_1perc[,2:3],1,mean)

#SHI01_01_003
SHI01_01_003 <- shingle_1perc[c(1,4:5)]
names(SHI01_01_003)[2] <- "PCR.1"
names(SHI01_01_003)[3] <- "PCR.2"

#add column of fish ID
SHI01_01_003$Sample <- "SHI01_01_003"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
SHI01_01_003_1perc <- subset(SHI01_01_003, PCR.1>=(sum(SHI01_01_003$PCR.1)*0.01) | PCR.2>=(sum(SHI01_01_003$PCR.2)*0.01))

#add empty PCR.3 so can bind all dfs together later
SHI01_01_003_1perc$PCR.3 <- NA

#reorganize to match other fish dfs
SHI01_01_003_1perc <- SHI01_01_003_1perc %>% relocate(PCR.3, .after = PCR.2)

#add total depth and average depth columns
SHI01_01_003_1perc$Total.depth <-apply(SHI01_01_003_1perc[,2:3],1,sum)
SHI01_01_003_1perc$Average.depth <-apply(SHI01_01_003_1perc[,2:3],1,mean)

#SHI02_01_001
SHI02_01_001 <- shingle_1perc[c(1,6:8)]
names(SHI02_01_001)[2] <- "PCR.1"
names(SHI02_01_001)[3] <- "PCR.2"
names(SHI02_01_001)[4] <- "PCR.3"

#add column of fish ID
SHI02_01_001$Sample <- "SHI02_01_001"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
SHI02_01_001_1perc <- subset(SHI02_01_001, PCR.1>=(sum(SHI02_01_001$PCR.1)*0.01) | PCR.2>=(sum(SHI02_01_001$PCR.2)*0.01) | PCR.3>=(sum(SHI02_01_001$PCR.3)*0.01))

#add total depth and average depth columns
SHI02_01_001_1perc$Total.depth <-apply(SHI02_01_001_1perc[,2:4],1,sum)
SHI02_01_001_1perc$Average.depth <-apply(SHI02_01_001_1perc[,2:4],1,mean)

#SHI02_01_002
SHI02_01_002 <- shingle_1perc[c(1,9:11)]
names(SHI02_01_002)[2] <- "PCR.1"
names(SHI02_01_002)[3] <- "PCR.2"
names(SHI02_01_002)[4] <- "PCR.3"

#add column of fish ID
SHI02_01_002$Sample <- "SHI02_01_002"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
SHI02_01_002_1perc <- subset(SHI02_01_002, PCR.1>=(sum(SHI02_01_002$PCR.1)*0.01) | PCR.2>=(sum(SHI02_01_002$PCR.2)*0.01) | PCR.3>=(sum(SHI02_01_002$PCR.3)*0.01))

#add total depth and average depth columns
SHI02_01_002_1perc$Total.depth <-apply(SHI02_01_002_1perc[,2:4],1,sum)
SHI02_01_002_1perc$Average.depth <-apply(SHI02_01_002_1perc[,2:4],1,mean)

#SHI02_01_003
SHI02_01_003 <- shingle_1perc[c(1,12:14)]
names(SHI02_01_003)[2] <- "PCR.1"
names(SHI02_01_003)[3] <- "PCR.2"
names(SHI02_01_003)[4] <- "PCR.3"

#add column of fish ID
SHI02_01_003$Sample <- "SHI02_01_003"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
SHI02_01_003_1perc <- subset(SHI02_01_003, PCR.1>=(sum(SHI02_01_003$PCR.1)*0.01) | PCR.2>=(sum(SHI02_01_003$PCR.2)*0.01) | PCR.3>=(sum(SHI02_01_003$PCR.3)*0.01))

#add total depth and average depth columns
SHI02_01_003_1perc$Total.depth <-apply(SHI02_01_003_1perc[,2:4],1,sum)
SHI02_01_003_1perc$Average.depth <-apply(SHI02_01_003_1perc[,2:4],1,mean)



#put it all together
shingle_1perc <- rbind(SHI01_01_002_1perc, SHI01_01_003_1perc, 
                   SHI02_01_001_1perc, SHI02_01_002_1perc, SHI02_01_003_1perc)

length(unique(shingle_1perc$ASV))

#keep only ASVs >=1% total depth of retained ASVs (to remove the very low depth ones)
shi_depths <- aggregate(Total.depth ~ ASV, data=shingle_1perc, FUN=sum)
shi_total <- sum(shi_depths[2])

shi_depths$Keep <- ifelse(shi_depths$Total.depth >= (shi_total*0.01), "yes","no")

#make df of only ASVs to keep
shi_keep <- subset(shi_depths, Keep=="yes", c(ASV))

#keep only those ASVs
shingle_1perc_final <- subset(shingle_1perc, ASV %in% shi_keep$ASV)


#add sequences - make it new from shingle folder
library (devtools)
library(tidyverse)
source_url("https://raw.githubusercontent.com/lrjoshi/FastaTabular/master/fasta_and_tabular.R")

FastaToTabular("eDNA/shingle/output/ASV/derep.fasta")
#The output will be stored as dna_table.csv in the current directory ("Bioinformatics")

#read in data
UNIXseqs <- read.csv("dna_table.csv")

#reformat, remove first column, name new 1st column to ASV and new 2nd to Sequence
#need the columns to match so can join to all_1perc 
UNIXseqs <- UNIXseqs[-1]
names(UNIXseqs)[1] <- "ASV"
names(UNIXseqs)[2] <- "Sequence"

#rename the contents of the ASV column to match the ASV column of all_1perc
#";.*" says: replace semicolon (;) and every character after that (.*), with nothing "".
UNIXseqs$ASV <- gsub(";.*", "", UNIXseqs$ASV)
#">" says: replace arrow (>) with nothing "".
UNIXseqs$ASV <- gsub(">", "", UNIXseqs$ASV)

#add data frames together
#library(tidyverse)
shingle_1perc_seqs <- shingle_1perc_final %>% 
  inner_join(UNIXseqs, by = c("ASV" = "ASV"))

#move sequence to the front
shingle_1perc_seqs <- shingle_1perc_seqs %>% relocate(Sequence, .after = ASV)

#add column of # PCR reps the reads came from
shingle_1perc_seqs$PCR.reps <- ifelse(shingle_1perc_seqs$PCR.1>0 & !is.na(shingle_1perc_seqs$PCR.3), 3, 2)


#add seq length
shingle_1perc_seqs$Sequence_length <- str_length(shingle_1perc_seqs$Sequence)

#reformat
shingle_1perc_seqs <- shingle_1perc_seqs %>% relocate(Sequence_length, .after = Sequence)
shi_depths_final <- aggregate(Total.depth ~ ASV, data=shingle_1perc_seqs, FUN=sum)


#save output for later analysis
write_csv(shingle_1perc_seqs, "eDNA/shingle/230510/shingleonly1%_ASVs.csv")




#release_02 (pos con removed)----
release <- all[-c(1:3)]

total_release <- sum(release[1:11])

#keep only samples that are >=1% of total release depth
release_1perc <- release[colSums(release[1:11]) >= (total_release*0.01)]

#add back ASV column
release_1perc <- cbind(ASV, release_1perc)


#pull out individual samples
#AWR000_02_001
AWR000_02_001 <- release_1perc[c(1:4)]
names(AWR000_02_001)[2] <- "PCR.1"
names(AWR000_02_001)[3] <- "PCR.2"
names(AWR000_02_001)[4] <- "PCR.3"

#add column of fish ID
AWR000_02_001$Sample <- "AWR000_02_001"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
AWR000_02_001_1perc <- subset(AWR000_02_001, PCR.1>=(sum(AWR000_02_001$PCR.1)*0.01) | PCR.2>=(sum(AWR000_02_001$PCR.2)*0.01) | PCR.3>=(sum(AWR000_02_001$PCR.3)*0.01))

#add total depth and average depth columns
AWR000_02_001_1perc$Total.depth <-apply(AWR000_02_001_1perc[,2:4],1,sum)
AWR000_02_001_1perc$Average.depth <-apply(AWR000_02_001_1perc[,2:4],1,mean)


#AWR000_02_002
AWR000_02_002 <- release_1perc[c(1, 5:7)]
names(AWR000_02_002)[2] <- "PCR.1"
names(AWR000_02_002)[3] <- "PCR.2"
names(AWR000_02_002)[4] <- "PCR.3"

#add column of fish ID
AWR000_02_002$Sample <- "AWR000_02_002"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
AWR000_02_002_1perc <- subset(AWR000_02_002, PCR.1>=(sum(AWR000_02_002$PCR.1)*0.01) | PCR.2>=(sum(AWR000_02_002$PCR.2)*0.01) | PCR.3>=(sum(AWR000_02_002$PCR.3)*0.01))

#add total depth and average depth columns
AWR000_02_002_1perc$Total.depth <-apply(AWR000_02_002_1perc[,2:4],1,sum)
AWR000_02_002_1perc$Average.depth <-apply(AWR000_02_002_1perc[,2:4],1,mean)


#AWR000_02_003
AWR000_02_003 <- release_1perc[c(1,8,9)]
names(AWR000_02_003)[2] <- "PCR.1"
names(AWR000_02_003)[3] <- "PCR.2"

#add column of fish ID
AWR000_02_003$Sample <- "AWR000_02_003"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
AWR000_02_003_1perc <- subset(AWR000_02_003, PCR.1>=(sum(AWR000_02_003$PCR.1)*0.01) | PCR.2>=(sum(AWR000_02_003$PCR.2)*0.01))

#add empty PCR.3 so can bind all dfs together later
AWR000_02_003_1perc$PCR.3 <- NA

#reorganize to match other fish dfs
AWR000_02_003_1perc <- AWR000_02_003_1perc %>% relocate(PCR.3, .after = PCR.2)

#add total depth and average depth columns
AWR000_02_003_1perc$Total.depth <-apply(AWR000_02_003_1perc[,2:3],1,sum)
AWR000_02_003_1perc$Average.depth <-apply(AWR000_02_003_1perc[,2:3],1,mean)



#put it all together
release_1perc <- rbind(AWR000_02_001_1perc, AWR000_02_002_1perc, AWR000_02_003_1perc)
length(unique(release_1perc$ASV))


#keep only ASVs >=1% total depth of retained ASVs (to remove the very low depth ones)
release_depths <- aggregate(Total.depth ~ ASV, data=release_1perc, FUN=sum)
release_total <- sum(release_depths[2])

release_depths$Keep <- ifelse(release_depths$Total.depth >= (release_total*0.01), "yes","no")

#make df of only ASVs to keep
release_keep <- subset(release_depths, Keep=="yes", c(ASV))

#keep only those ASVs
release_1perc_final <- subset(release_1perc, ASV %in% release_keep$ASV)
length(unique(release_1perc_final$ASV))
#9

#add sequences - make it new from release folder
library (devtools)
library(tidyverse)
source_url("https://raw.githubusercontent.com/lrjoshi/FastaTabular/master/fasta_and_tabular.R")

FastaToTabular("eDNA/release_02/output/ASV/derep.fasta")
#The output will be stored as dna_table.csv in the current directory ("UNIX-opt-new")
#rename the file so its clear it came from the release samples
file.rename("dna_table.csv", "dna_table_release_02.csv")

#read in data
UNIXseqs <- read.csv("dna_table_release_02.csv")


#reformat, remove first column, name new 1st column to ASV and new 2nd to Sequence
#need the columns to match so can join to all_1perc 
UNIXseqs <- UNIXseqs[-1]
names(UNIXseqs)[1] <- "ASV"
names(UNIXseqs)[2] <- "Sequence"

#rename the contents of the ASV column to match the ASV column of all_1perc
#";.*" says: replace semicolon (;) and every character after that (.*), with nothing "".
UNIXseqs$ASV <- gsub(";.*", "", UNIXseqs$ASV)
#">" says: replace arrow (>) with nothing "".
UNIXseqs$ASV <- gsub(">", "", UNIXseqs$ASV)

#add data frames together
#library(tidyverse)
release_1perc_seqs <- release_1perc_final %>% 
  inner_join(UNIXseqs, by = c("ASV" = "ASV"))

#move sequence to the front
release_1perc_seqs <- release_1perc_seqs %>% relocate(Sequence, .after = ASV)

#add column of # PCR reps the reads came from
release_1perc_seqs$PCR.reps <- ifelse(release_1perc_seqs$PCR.1>0 & !is.na(release_1perc_seqs$PCR.3), 3, 2)


#add seq length
release_1perc_seqs$Sequence_length <- str_length(release_1perc_seqs$Sequence)

#reformat
release_1perc_seqs <- release_1perc_seqs %>% relocate(Sequence_length, .after = Sequence)
release_depths_final <- aggregate(Total.depth ~ ASV, data=release_1perc_seqs, FUN=sum)


#save output for later analysis
write_csv(release_1perc_seqs, "eDNA/release_02/230512/releaseonly1%_ASVs.csv")




