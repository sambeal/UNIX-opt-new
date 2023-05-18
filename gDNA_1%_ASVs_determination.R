#analysis of the ASVs returned from the UNIX pipeline when run with default parameters except:
#sequence length restricted to 70-105bp (as determined from multiqc report)
#phred score 26 enforced during trimming
#sequences only taken through the dereplication stage

#230508 1% ASV determination
#assess holistically
#only keep samples (individual PCR reps) that >=1% total retained read depth
#to remove low quality PCRs while maintaining high quality ones

#managed to lose 230508 data, reassessed on 230516

#set working directory----
getwd()
#"/Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-opt-new"

#gDNA (230516 - e0.1)----
#load data
library(readr)
all <- readr::read_tsv("gDNA/output/ASV/derep.tsv", show_col_types = FALSE)
all <- as.data.frame(all)


#rename cols to be only fish name + index combo
names(all) = gsub(pattern = "_.*", replacement = "", x = names(all))
names(all)[1] <- "ASV"

#it's being weird when the ASV column is still there so will remove it and add it back before subsetting/fish
ASV <- all[c(1)]
all <- all[-c(1)]

#find total depth across all samples
total <- sum(all[1:54])

#cutoff=
total*0.01

#keep only samples that are >=1% of total depth
all_1perc <- all[colSums(all[1:54]) >= (total*0.01)]

#add back ASV column
all_1perc <- cbind(ASV, all_1perc)


#pull out individual fish
#190222-A1
A1190222 <- all_1perc[c(1:4)]
names(A1190222)[2] <- "PCR.1"
names(A1190222)[3] <- "PCR.2"
names(A1190222)[4] <- "PCR.3"

#add column of fish ID
A1190222$Fish <- "190222-A1"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
A1190222_1perc <- subset(A1190222, PCR.1>=(sum(A1190222$PCR.1)*0.01) | PCR.2>=(sum(A1190222$PCR.2)*0.01) | PCR.3>=(sum(A1190222$PCR.3)*0.01))

#add total depth and average depth columns
A1190222_1perc$Total.depth <-apply(A1190222_1perc[,2:4],1,sum)
A1190222_1perc$Average.depth <-apply(A1190222_1perc[,2:4],1,mean)

#190222-A2
A2190222 <- all_1perc[c(1,5:7)]
names(A2190222)[2] <- "PCR.1"
names(A2190222)[3] <- "PCR.2"
names(A2190222)[4] <- "PCR.3"

#add column of fish ID
A2190222$Fish <- "190922-A2"

#make new df of the 1% ASVs in this fish
A2190222_1perc <- subset(A2190222, PCR.1>=(sum(A2190222$PCR.1)*0.01) | PCR.2>=(sum(A2190222$PCR.2)*0.01) | PCR.3>=(sum(A2190222$PCR.3)*0.01))

#add total depth and average depth columns
A2190222_1perc$Total.depth <-apply(A2190222_1perc[,2:4],1,sum)
A2190222_1perc$Average.depth <-apply(A2190222_1perc[,2:4],1,mean)

#190222-M1
M1190222 <- all_1perc[c(1,8:10)]
names(M1190222)[2] <- "PCR.1"
names(M1190222)[3] <- "PCR.2"
names(M1190222)[4] <- "PCR.3"

#add column of fish ID
M1190222$Fish <- "190922-M1"

#make new df of the 1% ASVs in this fish
M1190222_1perc <- subset(M1190222, PCR.1>=(sum(M1190222$PCR.1)*0.01) | PCR.2>=(sum(M1190222$PCR.2)*0.01) | PCR.3>=(sum(M1190222$PCR.3)*0.01))

#add total depth and average depth columns
M1190222_1perc$Total.depth <-apply(M1190222_1perc[,2:4],1,sum)
M1190222_1perc$Average.depth <-apply(M1190222_1perc[,2:4],1,mean)

#190918-A1
A1190918 <- all_1perc[c(1,11:12)]
names(A1190918)[2] <- "PCR.1"
names(A1190918)[3] <- "PCR.2"

#add column of fish ID
A1190918$Fish <- "190918-A1"

#make new df of the 1% ASVs in this fish
A1190918_1perc <- subset(A1190918, PCR.1>=(sum(A1190918$PCR.1)*0.01) | PCR.2>=(sum(A1190918$PCR.2)*0.01))

#add PCR.3 column of NA's so can rbind all the dataframes together (doesnt work with diff col numbers)
A1190918_1perc$PCR.3 <- NA

#reorganize to match other fish dfs
library(dplyr)
A1190918_1perc <- A1190918_1perc %>% relocate(PCR.3, .after = PCR.2)

#add total depth and average depth columns
A1190918_1perc$Total.depth <-apply(A1190918_1perc[,2:3],1,sum)
A1190918_1perc$Average.depth <-apply(A1190918_1perc[,2:3],1,mean)


#190918-A3
A3190918 <- all_1perc[c(1,13:15)]
names(A3190918)[2] <- "PCR.1"
names(A3190918)[3] <- "PCR.2"
names(A3190918)[4] <- "PCR.3"

#add column of fish ID
A3190918$Fish <- "190918-A3"

#make new df of the 1% ASVs in this fish
A3190918_1perc <- subset(A3190918, PCR.1>=(sum(A3190918$PCR.1)*0.01) | PCR.2>=(sum(A3190918$PCR.2)*0.01) | PCR.3>=(sum(A3190918$PCR.3)*0.01))

#add total depth and average depth columns
A3190918_1perc$Total.depth <-apply(A3190918_1perc[,2:4],1,sum)
A3190918_1perc$Average.depth <-apply(A3190918_1perc[,2:4],1,mean)

#190918-A4
A4190918 <- all_1perc[c(1,16:18)]
names(A4190918)[2] <- "PCR.1"
names(A4190918)[3] <- "PCR.2"
names(A4190918)[4] <- "PCR.3"

#add column of fish ID
A4190918$Fish <- "190918-A4"

#make new df of the 1% ASVs in this fish
A4190918_1perc <- subset(A4190918, PCR.1>=(sum(A4190918$PCR.1)*0.01) | PCR.2>=(sum(A4190918$PCR.2)*0.01) | PCR.3>=(sum(A4190918$PCR.3)*0.01))

#add total depth and average depth columns
A4190918_1perc$Total.depth <-apply(A4190918_1perc[,2:4],1,sum)
A4190918_1perc$Average.depth <-apply(A4190918_1perc[,2:4],1,mean)

#190918-C4
C4190918 <- all_1perc[c(1,19:21)]
names(C4190918)[2] <- "PCR.1"
names(C4190918)[3] <- "PCR.2"
names(C4190918)[4] <- "PCR.3"

#add column of fish ID
C4190918$Fish <- "190918-C4"

#make new df of the 1% ASVs in this fish
C4190918_1perc <- subset(C4190918, PCR.1>=(sum(C4190918$PCR.1)*0.01) | PCR.2>=(sum(C4190918$PCR.2)*0.01) | PCR.3>=(sum(C4190918$PCR.3)*0.01))

#add total depth and average depth columns
C4190918_1perc$Total.depth <-apply(C4190918_1perc[,2:4],1,sum)
C4190918_1perc$Average.depth <-apply(C4190918_1perc[,2:4],1,mean)

#190918-D4
D4190918 <- all_1perc[c(1,22:24)]
names(D4190918)[2] <- "PCR.1"
names(D4190918)[3] <- "PCR.2"
names(D4190918)[4] <- "PCR.3"

#add column of fish ID
D4190918$Fish <- "190918-D4"

#make new df of the 1% ASVs in this fish
D4190918_1perc <- subset(D4190918, PCR.1>=(sum(D4190918$PCR.1)*0.01) | PCR.2>=(sum(D4190918$PCR.2)*0.01) | PCR.3>=(sum(D4190918$PCR.3)*0.01))

#add total depth and average depth columns
D4190918_1perc$Total.depth <-apply(D4190918_1perc[,2:4],1,sum)
D4190918_1perc$Average.depth <-apply(D4190918_1perc[,2:4],1,mean)

#190918-E7
E7190918 <- all_1perc[c(1,25:27)]
names(E7190918)[2] <- "PCR.1"
names(E7190918)[3] <- "PCR.2"
names(E7190918)[4] <- "PCR.3"

#add column of fish ID
E7190918$Fish <- "190918-E7"

#make new df of the 1% ASVs in this fish
E7190918_1perc <- subset(E7190918, PCR.1>=(sum(E7190918$PCR.1)*0.01) | PCR.2>=(sum(E7190918$PCR.2)*0.01) | PCR.3>=(sum(E7190918$PCR.3)*0.01))

#add total depth and average depth columns
E7190918_1perc$Total.depth <-apply(E7190918_1perc[,2:4],1,sum)
E7190918_1perc$Average.depth <-apply(E7190918_1perc[,2:4],1,mean)

#190918-F2
F2190918 <- all_1perc[c(1,28:30)]
names(F2190918)[2] <- "PCR.1"
names(F2190918)[3] <- "PCR.2"
names(F2190918)[4] <- "PCR.3"

#add column of fish ID
F2190918$Fish <- "190918-F2"

#make new df of the 1% ASVs in this fish
F2190918_1perc <- subset(F2190918, PCR.1>=(sum(F2190918$PCR.1)*0.01) | PCR.2>=(sum(F2190918$PCR.2)*0.01) | PCR.3>=(sum(F2190918$PCR.3)*0.01))

#add total depth and average depth columns
F2190918_1perc$Total.depth <-apply(F2190918_1perc[,2:4],1,sum)
F2190918_1perc$Average.depth <-apply(F2190918_1perc[,2:4],1,mean)

#190922-B2
B2190922 <- all_1perc[c(1,31:32)]
names(B2190922)[2] <- "PCR.1"
names(B2190922)[3] <- "PCR.2"

#add column of fish ID
B2190922$Fish <- "190922-B2"

#make new df of the 1% ASVs in this fish
B2190922_1perc <- subset(B2190922, PCR.1>=(sum(B2190922$PCR.1)*0.01) | PCR.2>=(sum(B2190922$PCR.2)*0.01))

#add PCR.3 column of NA's so can rbind all the dataframes together (doesnt work with diff col numbers)
B2190922_1perc$PCR.3 <- NA

#reorganize to match other fish dfs
library(dplyr)
B2190922_1perc <- B2190922_1perc %>% relocate(PCR.3, .after = PCR.2)

#add total depth and average depth columns
B2190922_1perc$Total.depth <-apply(B2190922_1perc[,2:3],1,sum)
B2190922_1perc$Average.depth <-apply(B2190922_1perc[,2:3],1,mean)

#190922-D2
D2190922 <- all_1perc[c(1,33:35)]
names(D2190922)[2] <- "PCR.1"
names(D2190922)[3] <- "PCR.2"
names(D2190922)[4] <- "PCR.3"

#add column of fish ID
D2190922$Fish <- "190922-D2"

#make new df of the 1% ASVs in this fish
D2190922_1perc <- subset(D2190922, PCR.1>=(sum(D2190922$PCR.1)*0.01) | PCR.2>=(sum(D2190922$PCR.2)*0.01) | PCR.3>=(sum(D2190922$PCR.3)*0.01))

#add total depth and average depth columns
D2190922_1perc$Total.depth <-apply(D2190922_1perc[,2:4],1,sum)
D2190922_1perc$Average.depth <-apply(D2190922_1perc[,2:4],1,mean)

#190922-F2
F2190922 <- all_1perc[c(1,36:38)]
names(F2190922)[2] <- "PCR.1"
names(F2190922)[3] <- "PCR.2"
names(F2190922)[4] <- "PCR.3"

#add column of fish ID
F2190922$Fish <- "190922-F2"

#make new df of the 1% ASVs in this fish
F2190922_1perc <- subset(F2190922, PCR.1>=(sum(F2190922$PCR.1)*0.01) | PCR.2>=(sum(F2190922$PCR.2)*0.01) | PCR.3>=(sum(F2190922$PCR.3)*0.01))

#add total depth and average depth columns
F2190922_1perc$Total.depth <-apply(F2190922_1perc[,2:4],1,sum)
F2190922_1perc$Average.depth <-apply(F2190922_1perc[,2:4],1,mean)

#190922-G2
G2190922 <- all_1perc[c(1,39:41)]
names(G2190922)[2] <- "PCR.1"
names(G2190922)[3] <- "PCR.2"
names(G2190922)[4] <- "PCR.3"

#add column of fish ID
G2190922$Fish <- "190922-G2"

#make new df of the 1% ASVs in this fish
G2190922_1perc <- subset(G2190922, PCR.1>=(sum(G2190922$PCR.1)*0.01) | PCR.2>=(sum(G2190922$PCR.2)*0.01) | PCR.3>=(sum(G2190922$PCR.3)*0.01))

#add total depth and average depth columns
G2190922_1perc$Total.depth <-apply(G2190922_1perc[,2:4],1,sum)
G2190922_1perc$Average.depth <-apply(G2190922_1perc[,2:4],1,mean)

#190922-H2
H2190922 <- all_1perc[c(1,42:44)]
names(H2190922)[2] <- "PCR.1"
names(H2190922)[3] <- "PCR.2"
names(H2190922)[4] <- "PCR.3"

#add column of fish ID
H2190922$Fish <- "190922-H2"

#make new df of the 1% ASVs in this fish
H2190922_1perc <- subset(H2190922, PCR.1>=(sum(H2190922$PCR.1)*0.01) | PCR.2>=(sum(H2190922$PCR.2)*0.01) | PCR.3>=(sum(H2190922$PCR.3)*0.01))

#add total depth and average depth columns
H2190922_1perc$Total.depth <-apply(H2190922_1perc[,2:4],1,sum)
H2190922_1perc$Average.depth <-apply(H2190922_1perc[,2:4],1,mean)

#190922-J2
J2190922 <- all_1perc[c(1,45:47)]
names(J2190922)[2] <- "PCR.1"
names(J2190922)[3] <- "PCR.2"
names(J2190922)[4] <- "PCR.3"

#add column of fish ID
J2190922$Fish <- "190922-J2"

#make new df of the 1% ASVs in this fish
J2190922_1perc <- subset(J2190922, PCR.1>=(sum(J2190922$PCR.1)*0.01) | PCR.2>=(sum(J2190922$PCR.2)*0.01) | PCR.3>=(sum(J2190922$PCR.3)*0.01))

#add total depth and average depth columns
J2190922_1perc$Total.depth <-apply(J2190922_1perc[,2:4],1,sum)
J2190922_1perc$Average.depth <-apply(J2190922_1perc[,2:4],1,mean)

#put it all together
all_1perc <- rbind(A1190222_1perc, A2190222_1perc, M1190222_1perc,
                   A1190918_1perc, A3190918_1perc, A4190918_1perc, C4190918_1perc,
                   D4190918_1perc, E7190918_1perc,F2190918_1perc,
                   B2190922_1perc, D2190922_1perc, F2190922_1perc,G2190922_1perc,H2190922_1perc,
                   J2190922_1perc)





#add sequences - make it new from gDNA folder
library (devtools)
library(tidyverse)
source_url("https://raw.githubusercontent.com/lrjoshi/FastaTabular/master/fasta_and_tabular.R")

FastaToTabular("gDNA/output/ASV/derep.fasta")
#The output will be stored as dna_table.csv in the current directory ("Bioinformatics")
file.rename("dna_table.csv", "dna_table_gDNA.csv")

#read in data
UNIXseqs <- read.csv("dna_table_gDNA.csv")
#read in data (from previously made .csv of UNIX seqs - moved into current working directory)

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
library(tidyverse)
all_1perc_seqs <- all_1perc %>% 
  inner_join(UNIXseqs, by = c("ASV" = "ASV"))

#move sequence to the front
all_1perc_seqs <- all_1perc_seqs %>% relocate(Sequence, .after = ASV)

#add column of # PCR reps the reads came from
all_1perc_seqs$PCR.reps <- ifelse(all_1perc_seqs$PCR.1>0 & !is.na(all_1perc_seqs$PCR.3), 3, 2)


#add seq length
library(tidyverse)
all_1perc_seqs$Sequence_length <- str_length(all_1perc_seqs$Sequence)

#reformat
all_1perc_seqs <- all_1perc_seqs %>% relocate(Sequence_length, .after = Sequence)


#save output for later analysis
#make folder with today's date
write_csv(all_1perc_seqs, "gDNA/230516/1%_ASVs.csv")




#gDNA (230518 - e0.2)----
#set working directory
getwd()
#"/Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-opt-new"

#load data
library(readr)
all <- readr::read_tsv("gDNA/output-e02/ASV/derep.tsv", show_col_types = FALSE)
all <- as.data.frame(all)


#rename cols to be only fish name + index combo
names(all) = gsub(pattern = "_.*", replacement = "", x = names(all))
names(all)[1] <- "ASV"

#it's being weird when the ASV column is still there so will remove it and add it back before subsetting/fish
ASV <- all[c(1)]
all <- all[-c(1)]

#find total depth across all samples
total <- sum(all[1:54])

#cutoff=
total*0.01

#keep only samples that are >=1% of total depth
all_1perc <- all[colSums(all[1:54]) >= (total*0.01)]

#add back ASV column
all_1perc <- cbind(ASV, all_1perc)


#pull out individual fish
#190222-A1
A1190222 <- all_1perc[c(1:4)]
names(A1190222)[2] <- "PCR.1"
names(A1190222)[3] <- "PCR.2"
names(A1190222)[4] <- "PCR.3"

#add column of fish ID
A1190222$Fish <- "190222-A1"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
A1190222_1perc <- subset(A1190222, PCR.1>=(sum(A1190222$PCR.1)*0.01) | PCR.2>=(sum(A1190222$PCR.2)*0.01) | PCR.3>=(sum(A1190222$PCR.3)*0.01))

#add total depth and average depth columns
A1190222_1perc$Total.depth <-apply(A1190222_1perc[,2:4],1,sum)
A1190222_1perc$Average.depth <-apply(A1190222_1perc[,2:4],1,mean)

#190222-A2
A2190222 <- all_1perc[c(1,5:7)]
names(A2190222)[2] <- "PCR.1"
names(A2190222)[3] <- "PCR.2"
names(A2190222)[4] <- "PCR.3"

#add column of fish ID
A2190222$Fish <- "190922-A2"

#make new df of the 1% ASVs in this fish
A2190222_1perc <- subset(A2190222, PCR.1>=(sum(A2190222$PCR.1)*0.01) | PCR.2>=(sum(A2190222$PCR.2)*0.01) | PCR.3>=(sum(A2190222$PCR.3)*0.01))

#add total depth and average depth columns
A2190222_1perc$Total.depth <-apply(A2190222_1perc[,2:4],1,sum)
A2190222_1perc$Average.depth <-apply(A2190222_1perc[,2:4],1,mean)

#190222-M1
M1190222 <- all_1perc[c(1,8:10)]
names(M1190222)[2] <- "PCR.1"
names(M1190222)[3] <- "PCR.2"
names(M1190222)[4] <- "PCR.3"

#add column of fish ID
M1190222$Fish <- "190922-M1"

#make new df of the 1% ASVs in this fish
M1190222_1perc <- subset(M1190222, PCR.1>=(sum(M1190222$PCR.1)*0.01) | PCR.2>=(sum(M1190222$PCR.2)*0.01) | PCR.3>=(sum(M1190222$PCR.3)*0.01))

#add total depth and average depth columns
M1190222_1perc$Total.depth <-apply(M1190222_1perc[,2:4],1,sum)
M1190222_1perc$Average.depth <-apply(M1190222_1perc[,2:4],1,mean)

#190918-A1
A1190918 <- all_1perc[c(1,11:12)]
names(A1190918)[2] <- "PCR.1"
names(A1190918)[3] <- "PCR.2"

#add column of fish ID
A1190918$Fish <- "190918-A1"

#make new df of the 1% ASVs in this fish
A1190918_1perc <- subset(A1190918, PCR.1>=(sum(A1190918$PCR.1)*0.01) | PCR.2>=(sum(A1190918$PCR.2)*0.01))

#add PCR.3 column of NA's so can rbind all the dataframes together (doesnt work with diff col numbers)
A1190918_1perc$PCR.3 <- NA

#reorganize to match other fish dfs
library(dplyr)
A1190918_1perc <- A1190918_1perc %>% relocate(PCR.3, .after = PCR.2)

#add total depth and average depth columns
A1190918_1perc$Total.depth <-apply(A1190918_1perc[,2:3],1,sum)
A1190918_1perc$Average.depth <-apply(A1190918_1perc[,2:3],1,mean)

#190918-A3
A3190918 <- all_1perc[c(1,13:15)]
names(A3190918)[2] <- "PCR.1"
names(A3190918)[3] <- "PCR.2"
names(A3190918)[4] <- "PCR.3"

#add column of fish ID
A3190918$Fish <- "190918-A3"

#make new df of the 1% ASVs in this fish
A3190918_1perc <- subset(A3190918, PCR.1>=(sum(A3190918$PCR.1)*0.01) | PCR.2>=(sum(A3190918$PCR.2)*0.01) | PCR.3>=(sum(A3190918$PCR.3)*0.01))

#add total depth and average depth columns
A3190918_1perc$Total.depth <-apply(A3190918_1perc[,2:4],1,sum)
A3190918_1perc$Average.depth <-apply(A3190918_1perc[,2:4],1,mean)

#190918-A4
A4190918 <- all_1perc[c(1,16:18)]
names(A4190918)[2] <- "PCR.1"
names(A4190918)[3] <- "PCR.2"
names(A4190918)[4] <- "PCR.3"

#add column of fish ID
A4190918$Fish <- "190918-A4"

#make new df of the 1% ASVs in this fish
A4190918_1perc <- subset(A4190918, PCR.1>=(sum(A4190918$PCR.1)*0.01) | PCR.2>=(sum(A4190918$PCR.2)*0.01) | PCR.3>=(sum(A4190918$PCR.3)*0.01))

#add total depth and average depth columns
A4190918_1perc$Total.depth <-apply(A4190918_1perc[,2:4],1,sum)
A4190918_1perc$Average.depth <-apply(A4190918_1perc[,2:4],1,mean)

#190918-C4
C4190918 <- all_1perc[c(1,19:21)]
names(C4190918)[2] <- "PCR.1"
names(C4190918)[3] <- "PCR.2"
names(C4190918)[4] <- "PCR.3"

#add column of fish ID
C4190918$Fish <- "190918-C4"

#make new df of the 1% ASVs in this fish
C4190918_1perc <- subset(C4190918, PCR.1>=(sum(C4190918$PCR.1)*0.01) | PCR.2>=(sum(C4190918$PCR.2)*0.01) | PCR.3>=(sum(C4190918$PCR.3)*0.01))

#add total depth and average depth columns
C4190918_1perc$Total.depth <-apply(C4190918_1perc[,2:4],1,sum)
C4190918_1perc$Average.depth <-apply(C4190918_1perc[,2:4],1,mean)

#190918-D4
D4190918 <- all_1perc[c(1,22:24)]
names(D4190918)[2] <- "PCR.1"
names(D4190918)[3] <- "PCR.2"
names(D4190918)[4] <- "PCR.3"

#add column of fish ID
D4190918$Fish <- "190918-D4"

#make new df of the 1% ASVs in this fish
D4190918_1perc <- subset(D4190918, PCR.1>=(sum(D4190918$PCR.1)*0.01) | PCR.2>=(sum(D4190918$PCR.2)*0.01) | PCR.3>=(sum(D4190918$PCR.3)*0.01))

#add total depth and average depth columns
D4190918_1perc$Total.depth <-apply(D4190918_1perc[,2:4],1,sum)
D4190918_1perc$Average.depth <-apply(D4190918_1perc[,2:4],1,mean)

#190918-E7
E7190918 <- all_1perc[c(1,25:27)]
names(E7190918)[2] <- "PCR.1"
names(E7190918)[3] <- "PCR.2"
names(E7190918)[4] <- "PCR.3"

#add column of fish ID
E7190918$Fish <- "190918-E7"

#make new df of the 1% ASVs in this fish
E7190918_1perc <- subset(E7190918, PCR.1>=(sum(E7190918$PCR.1)*0.01) | PCR.2>=(sum(E7190918$PCR.2)*0.01) | PCR.3>=(sum(E7190918$PCR.3)*0.01))

#add total depth and average depth columns
E7190918_1perc$Total.depth <-apply(E7190918_1perc[,2:4],1,sum)
E7190918_1perc$Average.depth <-apply(E7190918_1perc[,2:4],1,mean)

#190918-F2
F2190918 <- all_1perc[c(1,28:30)]
names(F2190918)[2] <- "PCR.1"
names(F2190918)[3] <- "PCR.2"
names(F2190918)[4] <- "PCR.3"

#add column of fish ID
F2190918$Fish <- "190918-F2"

#make new df of the 1% ASVs in this fish
F2190918_1perc <- subset(F2190918, PCR.1>=(sum(F2190918$PCR.1)*0.01) | PCR.2>=(sum(F2190918$PCR.2)*0.01) | PCR.3>=(sum(F2190918$PCR.3)*0.01))

#add total depth and average depth columns
F2190918_1perc$Total.depth <-apply(F2190918_1perc[,2:4],1,sum)
F2190918_1perc$Average.depth <-apply(F2190918_1perc[,2:4],1,mean)

#190922-B2
B2190922 <- all_1perc[c(1,31:32)]
names(B2190922)[2] <- "PCR.1"
names(B2190922)[3] <- "PCR.2"

#add column of fish ID
B2190922$Fish <- "190922-B2"

#make new df of the 1% ASVs in this fish
B2190922_1perc <- subset(B2190922, PCR.1>=(sum(B2190922$PCR.1)*0.01) | PCR.2>=(sum(B2190922$PCR.2)*0.01))

#add PCR.3 column of NA's so can rbind all the dataframes together (doesnt work with diff col numbers)
B2190922_1perc$PCR.3 <- NA

#reorganize to match other fish dfs
library(dplyr)
B2190922_1perc <- B2190922_1perc %>% relocate(PCR.3, .after = PCR.2)

#add total depth and average depth columns
B2190922_1perc$Total.depth <-apply(B2190922_1perc[,2:3],1,sum)
B2190922_1perc$Average.depth <-apply(B2190922_1perc[,2:3],1,mean)

#190922-D2
D2190922 <- all_1perc[c(1,33:35)]
names(D2190922)[2] <- "PCR.1"
names(D2190922)[3] <- "PCR.2"
names(D2190922)[4] <- "PCR.3"

#add column of fish ID
D2190922$Fish <- "190922-D2"

#make new df of the 1% ASVs in this fish
D2190922_1perc <- subset(D2190922, PCR.1>=(sum(D2190922$PCR.1)*0.01) | PCR.2>=(sum(D2190922$PCR.2)*0.01) | PCR.3>=(sum(D2190922$PCR.3)*0.01))

#add total depth and average depth columns
D2190922_1perc$Total.depth <-apply(D2190922_1perc[,2:4],1,sum)
D2190922_1perc$Average.depth <-apply(D2190922_1perc[,2:4],1,mean)

#190922-F2
F2190922 <- all_1perc[c(1,36:38)]
names(F2190922)[2] <- "PCR.1"
names(F2190922)[3] <- "PCR.2"
names(F2190922)[4] <- "PCR.3"

#add column of fish ID
F2190922$Fish <- "190922-F2"

#make new df of the 1% ASVs in this fish
F2190922_1perc <- subset(F2190922, PCR.1>=(sum(F2190922$PCR.1)*0.01) | PCR.2>=(sum(F2190922$PCR.2)*0.01) | PCR.3>=(sum(F2190922$PCR.3)*0.01))

#add total depth and average depth columns
F2190922_1perc$Total.depth <-apply(F2190922_1perc[,2:4],1,sum)
F2190922_1perc$Average.depth <-apply(F2190922_1perc[,2:4],1,mean)

#190922-G2
G2190922 <- all_1perc[c(1,39:41)]
names(G2190922)[2] <- "PCR.1"
names(G2190922)[3] <- "PCR.2"
names(G2190922)[4] <- "PCR.3"

#add column of fish ID
G2190922$Fish <- "190922-G2"

#make new df of the 1% ASVs in this fish
G2190922_1perc <- subset(G2190922, PCR.1>=(sum(G2190922$PCR.1)*0.01) | PCR.2>=(sum(G2190922$PCR.2)*0.01) | PCR.3>=(sum(G2190922$PCR.3)*0.01))

#add total depth and average depth columns
G2190922_1perc$Total.depth <-apply(G2190922_1perc[,2:4],1,sum)
G2190922_1perc$Average.depth <-apply(G2190922_1perc[,2:4],1,mean)

#190922-H2
H2190922 <- all_1perc[c(1,42:44)]
names(H2190922)[2] <- "PCR.1"
names(H2190922)[3] <- "PCR.2"
names(H2190922)[4] <- "PCR.3"

#add column of fish ID
H2190922$Fish <- "190922-H2"

#make new df of the 1% ASVs in this fish
H2190922_1perc <- subset(H2190922, PCR.1>=(sum(H2190922$PCR.1)*0.01) | PCR.2>=(sum(H2190922$PCR.2)*0.01) | PCR.3>=(sum(H2190922$PCR.3)*0.01))

#add total depth and average depth columns
H2190922_1perc$Total.depth <-apply(H2190922_1perc[,2:4],1,sum)
H2190922_1perc$Average.depth <-apply(H2190922_1perc[,2:4],1,mean)

#190922-J2
J2190922 <- all_1perc[c(1,45:47)]
names(J2190922)[2] <- "PCR.1"
names(J2190922)[3] <- "PCR.2"
names(J2190922)[4] <- "PCR.3"

#add column of fish ID
J2190922$Fish <- "190922-J2"

#make new df of the 1% ASVs in this fish
J2190922_1perc <- subset(J2190922, PCR.1>=(sum(J2190922$PCR.1)*0.01) | PCR.2>=(sum(J2190922$PCR.2)*0.01) | PCR.3>=(sum(J2190922$PCR.3)*0.01))

#add total depth and average depth columns
J2190922_1perc$Total.depth <-apply(J2190922_1perc[,2:4],1,sum)
J2190922_1perc$Average.depth <-apply(J2190922_1perc[,2:4],1,mean)

#put it all together
all_1perc <- rbind(A1190222_1perc, A2190222_1perc, M1190222_1perc,
                   A1190918_1perc, A3190918_1perc, A4190918_1perc, C4190918_1perc,
                   D4190918_1perc, E7190918_1perc,F2190918_1perc,
                   B2190922_1perc, D2190922_1perc, F2190922_1perc,G2190922_1perc,H2190922_1perc,
                   J2190922_1perc)





#add sequences - make it new from gDNA folder
library (devtools)
library(tidyverse)
source_url("https://raw.githubusercontent.com/lrjoshi/FastaTabular/master/fasta_and_tabular.R")

FastaToTabular("gDNA/output-e02/ASV/derep.fasta")
#The output will be stored as dna_table.csv in the current directory ("Bioinformatics")
file.rename("dna_table.csv", "dna_table_gDNAe02.csv")

#read in data
UNIXseqs <- read.csv("dna_table_gDNAe02.csv")
#read in data (from previously made .csv of UNIX seqs - moved into current working directory)

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
library(tidyverse)
all_1perc_seqs <- all_1perc %>% 
  inner_join(UNIXseqs, by = c("ASV" = "ASV"))

#move sequence to the front
all_1perc_seqs <- all_1perc_seqs %>% relocate(Sequence, .after = ASV)

#add column of # PCR reps the reads came from
all_1perc_seqs$PCR.reps <- ifelse(all_1perc_seqs$PCR.1>0 & !is.na(all_1perc_seqs$PCR.3), 3, 2)


#add seq length
library(tidyverse)
all_1perc_seqs$Sequence_length <- str_length(all_1perc_seqs$Sequence)

#reformat
all_1perc_seqs <- all_1perc_seqs %>% relocate(Sequence_length, .after = Sequence)


#save output for later analysis
#make folder with today's date
write_csv(all_1perc_seqs, "gDNA/230518/1%_ASVs.csv")



