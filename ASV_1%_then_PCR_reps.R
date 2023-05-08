#analysis of the ASVs returned from the UNIX pipeline when run with default parameters except:
#sequence length restricted to 70-105bp (as determined from multiqc report)
#phred score 26 enforced during trimming
#sequences only taken through the dereplication stage

#Parse through all ASVs to find only those with depths >=10% per fish

#set working directory----
getwd()

#load data----
library(readr)
all <-readr::read_tsv("output1/ASV/derep.tsv", show_col_types = FALSE)

#rename cols to be only fihs name
names(all) = gsub(pattern = "_.*", replacement = "", x = names(all))
names(all)[1] <- "ASV"

#pull out individual fish----
#190222-A1----
A1190222n1 <- all[c(1:4)]
names(A1190222n1)[1] <- "ASV"
names(A1190222n1)[2] <- "PCR.1"
names(A1190222n1)[3] <- "PCR.2"
names(A1190222n1)[4] <- "PCR.3"

#add column of fish ID
A1190222n1$Fish <- "190922-A1"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
A1190222_1perc <- subset(A1190222n1, PCR.1>=(sum(A1190222n1$PCR.1)*0.01) | PCR.2>=(sum(A1190222n1$PCR.2)*0.01) | PCR.3>=(sum(A1190222n1$PCR.3)*0.01))

#190222-A2----
A2190222n1 <- all[c(1,5:7)]
names(A2190222n1)[1] <- "ASV"
names(A2190222n1)[2] <- "PCR.1"
names(A2190222n1)[3] <- "PCR.2"
names(A2190222n1)[4] <- "PCR.3"

#add column of fish ID
A2190222n1$Fish <- "190922-A2"

#make new df of the 1% ASVs in this fish
A2190222_1perc <- subset(A2190222n1, PCR.1>=(sum(A2190222n1$PCR.1)*0.01) | PCR.2>=(sum(A2190222n1$PCR.2)*0.01) | PCR.3>=(sum(A2190222n1$PCR.3)*0.01))


#190222-H1----
H1190222n1 <- all[c(1,8:10)]
names(H1190222n1)[1] <- "ASV"
names(H1190222n1)[2] <- "PCR.1"
names(H1190222n1)[3] <- "PCR.2"
names(H1190222n1)[4] <- "PCR.3"

#add column of fish ID
H1190222n1$Fish <- "190922-H1"

#make new df of the 1% ASVs in this fish
H1190222_1perc <- subset(H1190222n1, PCR.1>=(sum(H1190222n1$PCR.1)*0.01) | PCR.2>=(sum(H1190222n1$PCR.2)*0.01) | PCR.3>=(sum(H1190222n1$PCR.3)*0.01))

#cut off for PCR.1 = 22.63, PCR.2=12.6,PCR.3=0.86
#i think PCR.3, even though removed later on, is why I'm having so many low-depth ASVs make it through to the end


#190222-M1----
M1190222n1 <- all[c(1,11:13)]
names(M1190222n1)[1] <- "ASV"
names(M1190222n1)[2] <- "PCR.1"
names(M1190222n1)[3] <- "PCR.2"
names(M1190222n1)[4] <- "PCR.3"

#add column of fish ID
M1190222n1$Fish <- "190922-M1"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
M1190222_1perc <- subset(M1190222n1, PCR.1>=(sum(M1190222n1$PCR.1)*0.01) | PCR.2>=(sum(M1190222n1$PCR.2)*0.01) | PCR.3>=(sum(M1190222n1$PCR.3)*0.01))

#190918-A1----
A1190918n1 <- all[c(1,14:16)]
names(A1190918n1)[1] <- "ASV"
names(A1190918n1)[2] <- "PCR.1"
names(A1190918n1)[3] <- "PCR.2"
names(A1190918n1)[4] <- "PCR.3"

#add column of fish ID
A1190918n1$Fish <- "190918-A1"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
A1190918_1perc <- subset(A1190918n1, PCR.1>=(sum(A1190918n1$PCR.1)*0.01) | PCR.2>=(sum(A1190918n1$PCR.2)*0.01) | PCR.3>=(sum(A1190918n1$PCR.3)*0.01))

#190918-A3----
A3190918n1 <- all[c(1,17:19)]
names(A3190918n1)[1] <- "ASV"
names(A3190918n1)[2] <- "PCR.1"
names(A3190918n1)[3] <- "PCR.2"
names(A3190918n1)[4] <- "PCR.3"

#add column of fish ID
A3190918n1$Fish <- "190918-A3"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
A3190918_1perc <- subset(A3190918n1, PCR.1>=(sum(A3190918n1$PCR.1)*0.01) | PCR.2>=(sum(A3190918n1$PCR.2)*0.01) | PCR.3>=(sum(A3190918n1$PCR.3)*0.01))

#190918-A4----
A4190918n1 <- all[c(1,20:22)]
names(A4190918n1)[1] <- "ASV"
names(A4190918n1)[2] <- "PCR.1"
names(A4190918n1)[3] <- "PCR.2"
names(A4190918n1)[4] <- "PCR.3"

#add column of fish ID
A4190918n1$Fish <- "190918-A4"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
A4190918_1perc <- subset(A4190918n1, PCR.1>=(sum(A4190918n1$PCR.1)*0.01) | PCR.2>=(sum(A4190918n1$PCR.2)*0.01) | PCR.3>=(sum(A4190918n1$PCR.3)*0.01))

#190918-C4----
C4190918n1 <- all[c(1,23:25)]
names(C4190918n1)[1] <- "ASV"
names(C4190918n1)[2] <- "PCR.1"
names(C4190918n1)[3] <- "PCR.2"
names(C4190918n1)[4] <- "PCR.3"

#add column of fish ID
C4190918n1$Fish <- "190918-C4"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
C4190918_1perc <- subset(C4190918n1, PCR.1>=(sum(C4190918n1$PCR.1)*0.01) | PCR.2>=(sum(C4190918n1$PCR.2)*0.01) | PCR.3>=(sum(C4190918n1$PCR.3)*0.01))

#190918-D4----
D4190918n1 <- all[c(1,26:28)]
names(D4190918n1)[1] <- "ASV"
names(D4190918n1)[2] <- "PCR.1"
names(D4190918n1)[3] <- "PCR.2"
names(D4190918n1)[4] <- "PCR.3"

#add column of fish ID
D4190918n1$Fish <- "190918-D4"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
D4190918_1perc <- subset(D4190918n1, PCR.1>=(sum(D4190918n1$PCR.1)*0.01) | PCR.2>=(sum(D4190918n1$PCR.2)*0.01) | PCR.3>=(sum(D4190918n1$PCR.3)*0.01))

#190918-E7----
E7190918n1 <- all[c(1,29:31)]
names(E7190918n1)[1] <- "ASV"
names(E7190918n1)[2] <- "PCR.1"
names(E7190918n1)[3] <- "PCR.2"
names(E7190918n1)[4] <- "PCR.3"

#add column of fish ID
E7190918n1$Fish <- "190918-E7"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
E7190918_1perc <- subset(E7190918n1, PCR.1>=(sum(E7190918n1$PCR.1)*0.01) | PCR.2>=(sum(E7190918n1$PCR.2)*0.01) | PCR.3>=(sum(E7190918n1$PCR.3)*0.01))

#190918-F2----
F2190918n1 <- all[c(1,32:34)]
names(F2190918n1)[1] <- "ASV"
names(F2190918n1)[2] <- "PCR.1"
names(F2190918n1)[3] <- "PCR.2"
names(F2190918n1)[4] <- "PCR.3"

#add column of fish ID
F2190918n1$Fish <- "190918-F2"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
F2190918_1perc <- subset(F2190918n1, PCR.1>=(sum(F2190918n1$PCR.1)*0.01) | PCR.2>=(sum(F2190918n1$PCR.2)*0.01) | PCR.3>=(sum(F2190918n1$PCR.3)*0.01))

#190922-B2----
B2190922n1 <- all[c(1,35:37)]
names(B2190922n1)[1] <- "ASV"
names(B2190922n1)[2] <- "PCR.1"
names(B2190922n1)[3] <- "PCR.2"
names(B2190922n1)[4] <- "PCR.3"

#add column of fish ID
B2190922n1$Fish <- "190922-B2"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
B2190922_1perc <- subset(B2190922n1, PCR.1>=(sum(B2190922n1$PCR.1)*0.01) | PCR.2>=(sum(B2190922n1$PCR.2)*0.01) | PCR.3>=(sum(B2190922n1$PCR.3)*0.01))

#190922-D2----
D2190922n1 <- all[c(1,38:40)]
names(D2190922n1)[1] <- "ASV"
names(D2190922n1)[2] <- "PCR.1"
names(D2190922n1)[3] <- "PCR.2"
names(D2190922n1)[4] <- "PCR.3"

#add column of fish ID
D2190922n1$Fish <- "190922-D2"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
D2190922_1perc <- subset(D2190922n1, PCR.1>=(sum(D2190922n1$PCR.1)*0.01) | PCR.2>=(sum(D2190922n1$PCR.2)*0.01) | PCR.3>=(sum(D2190922n1$PCR.3)*0.01))

#190922-F2----
F2190922n1 <- all[c(1,41:43)]
names(F2190922n1)[1] <- "ASV"
names(F2190922n1)[2] <- "PCR.1"
names(F2190922n1)[3] <- "PCR.2"
names(F2190922n1)[4] <- "PCR.3"

#add column of fish ID
F2190922n1$Fish <- "190922-F2"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
F2190922_1perc <- subset(F2190922n1, PCR.1>=(sum(F2190922n1$PCR.1)*0.01) | PCR.2>=(sum(F2190922n1$PCR.2)*0.01) | PCR.3>=(sum(F2190922n1$PCR.3)*0.01))

#190922-G2----
G2190922n1 <- all[c(1,44:46)]
names(G2190922n1)[1] <- "ASV"
names(G2190922n1)[2] <- "PCR.1"
names(G2190922n1)[3] <- "PCR.2"
names(G2190922n1)[4] <- "PCR.3"

#add column of fish ID
G2190922n1$Fish <- "190922-G2"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
G2190922_1perc <- subset(G2190922n1, PCR.1>=(sum(G2190922n1$PCR.1)*0.01) | PCR.2>=(sum(G2190922n1$PCR.2)*0.01) | PCR.3>=(sum(G2190922n1$PCR.3)*0.01))

#190922-H2----
H2190922n1 <- all[c(1,47:49)]
names(H2190922n1)[1] <- "ASV"
names(H2190922n1)[2] <- "PCR.1"
names(H2190922n1)[3] <- "PCR.2"
names(H2190922n1)[4] <- "PCR.3"

#add column of fish ID
H2190922n1$Fish <- "190922-H2"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
H2190922_1perc <- subset(H2190922n1, PCR.1>=(sum(H2190922n1$PCR.1)*0.01) | PCR.2>=(sum(H2190922n1$PCR.2)*0.01) | PCR.3>=(sum(H2190922n1$PCR.3)*0.01))

#190922-J2----
J2190922n1 <- all[c(1,50:52)]
names(J2190922n1)[1] <- "ASV"
names(J2190922n1)[2] <- "PCR.1"
names(J2190922n1)[3] <- "PCR.2"
names(J2190922n1)[4] <- "PCR.3"

#add column of fish ID
J2190922n1$Fish <- "190922-J2"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
J2190922_1perc <- subset(J2190922n1, PCR.1>=(sum(J2190922n1$PCR.1)*0.01) | PCR.2>=(sum(J2190922n1$PCR.2)*0.01) | PCR.3>=(sum(J2190922n1$PCR.3)*0.01))

#NTC----
NTCn1 <- all[c(1,53:55)]
names(NTCn1)[1] <- "ASV"
names(NTCn1)[2] <- "PCR.1"
names(NTCn1)[3] <- "PCR.2"
names(NTCn1)[4] <- "PCR.3"

#add column of fish ID
NTCn1$Fish <- "NTC"

#make new df of the 1% ASVs in this fish PER PCR REP (not across all reps)
NTC_1perc <- subset(NTCn1, PCR.1>=(sum(NTCn1$PCR.1)*0.01) | PCR.2>=(sum(NTCn1$PCR.2)*0.01) | PCR.3>=(sum(NTCn1$PCR.3)*0.01))





#put it all together----
all_1perc <- rbind(A1190222_1perc,A2190222_1perc, H1190222_1perc,M1190222_1perc,
                   A1190918_1perc,A3190918_1perc,A4190918_1perc,C4190918_1perc,
                   D4190918_1perc, E7190918_1perc,F2190918_1perc,B2190922_1perc,
                   D2190922_1perc, F2190922_1perc,G2190922_1perc,H2190922_1perc,
                   J2190922_1perc,NTC_1perc)

#add columns for total depth
all_1perc$Total.depth <-apply(all_1perc[,2:4],1,sum)

#rearrange
library(tidyverse)
all_1perc <- all_1perc %>% relocate(Total.depth, .after = PCR.3)


#add sequences
#read in data (from previously made .csv of UNIX seqs - moved into current working directory)
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
library(tidyverse)
all_1perc_seqs <- all_1perc %>% 
  inner_join(UNIXseqs, by = c("ASV" = "ASV"))

#move sequence to the front
all_1perc_seqs <- all_1perc_seqs %>% relocate(Sequence, .after = ASV)

#save output for later analysis
write_csv(all_1perc_seqs, "organized_ASVs/1%_ASVs.csv")




#assess if the ASV is also consistently detected across replicates----
all_1perc_seqs <- read.csv("organized_ASVs/1%_ASVs.csv")

#change data types
all_1perc_seqs$Fish <- as.factor(all_1perc_seqs$Fish)
all_1perc_seqs$ASV <- as.factor(all_1perc_seqs$ASV)
all_1perc_seqs$Sequence <- as.factor(all_1perc_seqs$Sequence)

#make a list of the fish, each in its own dataframe
fish_1perc <- split(all_1perc_seqs, all_1perc_seqs$Fish)
str(fish_1perc)
fish <- lapply(seq_along(fish_1perc), function(x) as.data.frame(fish_1perc[[x]])[,1:7])

#rename and make the data frames----
A1190918_1perc <- fish[[1]]
A3190918_1perc <- fish[[2]]
A4190918_1perc  <- fish[[3]]
C4190918_1perc <- fish[[4]]
D4190918_1perc  <- fish[[5]]
E7190918_1perc  <- fish[[6]]
F2190918_1perc  <- fish[[7]]
A1190222_1perc  <- fish[[8]]
A2190222_1perc <- fish[[9]]
B2190922_1perc  <- fish[[10]]
D2190922_1perc  <- fish[[11]]
F2190922_1perc  <- fish[[12]]
G2190922_1perc  <- fish[[13]]
H1190222_1perc  <- fish[[14]]
H2190922_1perc  <- fish[[15]]
J2190922_1perc <- fish[[16]]
M1190222_1perc <- fish[[17]]
NTC_1perc  <- fish[[18]]

#assess quality of PCR reps (ie. did all work?)----
#if the total depth of a PCR rep is less than 10% of the total depth of the fish, the rep will be removed entirely

#190918-A1----
#turn ASV into row names, remove non-numeric columns
A1190918_true <- A1190918_1perc[,-1]
rownames(A1190918_true) <- A1190918_1perc[,1]
A1190918_true <- A1190918_true[c(2:4)]

#subset to keep only PCR reps where the total depth > 10% of the reads
A1190918_true <- A1190918_true[,colSums(A1190918_true) >= (sum(A1190918_1perc$Total.depth)*0.1)]

#get new total depth of these ASVs
A1190918_true$Total.depth <- apply(A1190918_true[,1:2],1,sum)

#add column for number of PCR reps this total depth came from
A1190918_true$Num.PCR.reps <- "2"

#add back in fish name
A1190918_true$Fish <- "190918-A1"

#make ASV the first column again
A1190918_true <- tibble::rownames_to_column(A1190918_true, "ASV")

#save new dataframe of only ASV, total depth, and fish
A1190918_final <- A1190918_true[c(1,4:6)]

#190918-A3----
#turn ASV into row names, remove non-numeric columns
A3190918_true <- A3190918_1perc[,-1]
rownames(A3190918_true) <- A3190918_1perc[,1]
A3190918_true <- A3190918_true[c(2:4)]

#subset to keep only PCR reps where the total depth > 10% of the reads
A3190918_true <- A3190918_true[,colSums(A3190918_true) >= (sum(A3190918_1perc$Total.depth)*0.1)]

#get new total depth of these ASVs
A3190918_true$Total.depth <- apply(A3190918_true[,1:3],1,sum)

#add column for number of PCR reps this total depth came from
A3190918_true$Num.PCR.reps <- "3"

#add back in fish name
A3190918_true$Fish <- "190918-A3"

#make ASV the first column again
A3190918_true <- tibble::rownames_to_column(A3190918_true, "ASV")

#save new dataframe of only ASV, total depth, and fish
A3190918_final <- A3190918_true[c(1,5:7)]

#190918-A4----
#turn ASV into row names, remove non-numeric columns
A4190918_true <- A4190918_1perc[,-1]
rownames(A4190918_true) <- A4190918_1perc[,1]
A4190918_true <- A4190918_true[c(2:4)]

#subset to keep only PCR reps where the total depth > 10% of the reads
A4190918_true <- A4190918_true[,colSums(A4190918_true) >= (sum(A4190918_1perc$Total.depth)*0.1)]

#get new total depth of these ASVs
A4190918_true$Total.depth <- apply(A4190918_true[,1:3],1,sum)

#add column for number of PCR reps this total depth came from
A4190918_true$Num.PCR.reps <- "3"

#add back in fish name
A4190918_true$Fish <- "190918-A4"

#make ASV the first column again
A4190918_true <- tibble::rownames_to_column(A4190918_true, "ASV")

#save new dataframe of only ASV, total depth, and fish
A4190918_final <- A4190918_true[c(1,5:7)]

#190918-C4----
#turn ASV into row names, remove non-numeric columns
C4190918_true <- C4190918_1perc[,-1]
rownames(C4190918_true) <- C4190918_1perc[,1]
C4190918_true <- C4190918_true[c(2:4)]

#subset to keep only PCR reps where the total depth > 10% of the reads
C4190918_true <- C4190918_true[,colSums(C4190918_true) >= (sum(C4190918_1perc$Total.depth)*0.1)]

#get new total depth of these ASVs
C4190918_true$Total.depth <- apply(C4190918_true[,1:3],1,sum)

#add column for number of PCR reps this total depth came from
C4190918_true$Num.PCR.reps <- "3"

#add back in fish name
C4190918_true$Fish <- "190918-C4"

#make ASV the first column again
C4190918_true <- tibble::rownames_to_column(C4190918_true, "ASV")

#save new dataframe of only ASV, total depth, and fish
C4190918_final <- C4190918_true[c(1,5:7)]

#190918-D4----
#turn ASV into row names, remove non-numeric columns
D4190918_true <- D4190918_1perc[,-1]
rownames(D4190918_true) <- D4190918_1perc[,1]
D4190918_true <- D4190918_true[c(2:4)]

#subset to keep only PCR reps where the total depth > 10% of the reads
D4190918_true <- D4190918_true[,colSums(D4190918_true) >= (sum(D4190918_1perc$Total.depth)*0.1)]

#get new total depth of these ASVs
D4190918_true$Total.depth <- apply(D4190918_true[,1:3],1,sum)

#add column for number of PCR reps this total depth came from
D4190918_true$Num.PCR.reps <- "3"

#add back in fish name
D4190918_true$Fish <- "190918-D4"

#make ASV the first column again
D4190918_true <- tibble::rownames_to_column(D4190918_true, "ASV")

#save new dataframe of only ASV, total depth, and fish
D4190918_final <- D4190918_true[c(1,5:7)]

#190918-E7----
#turn ASV into row names, remove non-numeric columns
E7190918_true <- E7190918_1perc[,-1]
rownames(E7190918_true) <- E7190918_1perc[,1]
E7190918_true <- E7190918_true[c(2:4)]

#subset to keep only PCR reps where the total depth > 10% of the reads
E7190918_true <- E7190918_true[,colSums(E7190918_true) >= (sum(E7190918_1perc$Total.depth)*0.1)]

#get new total depth of these ASVs
E7190918_true$Total.depth <- apply(E7190918_true[,1:3],1,sum)

#add column for number of PCR reps this total depth came from
E7190918_true$Num.PCR.reps <- "3"

#add back in fish name
E7190918_true$Fish <- "190918-E7"

#make ASV the first column again
E7190918_true <- tibble::rownames_to_column(E7190918_true, "ASV")

#save new dataframe of only ASV, total depth, and fish
E7190918_final <- E7190918_true[c(1,5:7)]

#190918-F2----
#turn ASV into row names, remove non-numeric columns
F2190918_true <- F2190918_1perc[,-1]
rownames(F2190918_true) <- F2190918_1perc[,1]
F2190918_true <- F2190918_true[c(2:4)]

#subset to keep only PCR reps where the total depth > 10% of the reads
F2190918_true <- F2190918_true[,colSums(F2190918_true) >= (sum(F2190918_1perc$Total.depth)*0.1)]

#get new total depth of these ASVs
F2190918_true$Total.depth <- apply(F2190918_true[,1:3],1,sum)

#add column for number of PCR reps this total depth came from
F2190918_true$Num.PCR.reps <- "3"

#add back in fish name
F2190918_true$Fish <- "190918-F2"

#make ASV the first column again
F2190918_true <- tibble::rownames_to_column(F2190918_true, "ASV")

#save new dataframe of only ASV, total depth, and fish
F2190918_final <- F2190918_true[c(1,5:7)]

#190222-A1----
#turn ASV into row names, remove non-numeric columns
A1190222_true <- A1190222_1perc[,-1]
rownames(A1190222_true) <- A1190222_1perc[,1]
A1190222_true <- A1190222_true[c(2:4)]

#subset to keep only PCR reps where the total depth > 10% of the reads
A1190222_true <- A1190222_true[,colSums(A1190222_true) >= (sum(A1190222_1perc$Total.depth)*0.1)]

#get new total depth of these ASVs
A1190222_true$Total.depth <- apply(A1190222_true[,1:3],1,sum)

#add column for number of PCR reps this total depth came from
A1190222_true$Num.PCR.reps <- "3"

#add back in fish name
A1190222_true$Fish <- "190222-A1"

#make ASV the first column again
A1190222_true <- tibble::rownames_to_column(A1190222_true, "ASV")

#save new dataframe of only ASV, total depth, and fish
A1190222_final <- A1190222_true[c(1,5:7)]

#190222-A2----
#turn ASV into row names, remove non-numeric columns
A2190222_true <- A2190222_1perc[,-1]
rownames(A2190222_true) <- A2190222_1perc[,1]
A2190222_true <- A2190222_true[c(2:4)]

#subset to keep only PCR reps where the total depth > 10% of the reads
A2190222_true <- A2190222_true[,colSums(A2190222_true) >= (sum(A2190222_1perc$Total.depth)*0.1)]

#get new total depth of these ASVs
A2190222_true$Total.depth <- apply(A2190222_true[,1:3],1,sum)

#add column for number of PCR reps this total depth came from
A2190222_true$Num.PCR.reps <- "3"

#add back in fish name
A2190222_true$Fish <- "190222-A2"

#make ASV the first column again
A2190222_true <- tibble::rownames_to_column(A2190222_true, "ASV")

#save new dataframe of only ASV, total depth, and fish
A2190222_final <- A2190222_true[c(1,5:7)]

#190922-B2----
#turn ASV into row names, remove non-numeric columns
B2190922_true <- B2190922_1perc[,-1]
rownames(B2190922_true) <- B2190922_1perc[,1]
B2190922_true <- B2190922_true[c(2:4)]

#subset to keep only PCR reps where the total depth > 10% of the reads
B2190922_true <- B2190922_true[,colSums(B2190922_true) >= (sum(B2190922_1perc$Total.depth)*0.1)]

#get new total depth of these ASVs
B2190922_true$Total.depth <- apply(B2190922_true[,1:2],1,sum)

#add column for number of PCR reps this total depth came from
B2190922_true$Num.PCR.reps <- "2"

#add back in fish name
B2190922_true$Fish <- "190922-B2"

#make ASV the first column again
B2190922_true <- tibble::rownames_to_column(B2190922_true, "ASV")

#save new dataframe of only ASV, total depth, and fish
B2190922_final <- B2190922_true[c(1,4:6)]

#190922-D2----
#turn ASV into row names, remove non-numeric columns
D2190922_true <- D2190922_1perc[,-1]
rownames(D2190922_true) <- D2190922_1perc[,1]
D2190922_true <- D2190922_true[c(2:4)]

#subset to keep only PCR reps where the total depth > 10% of the reads
D2190922_true <- D2190922_true[,colSums(D2190922_true) >= (sum(D2190922_1perc$Total.depth)*0.1)]

#get new total depth of these ASVs
D2190922_true$Total.depth <- apply(D2190922_true[,1:3],1,sum)

#add column for number of PCR reps this total depth came from
D2190922_true$Num.PCR.reps <- "3"

#add back in fish name
D2190922_true$Fish <- "190922-D2"

#make ASV the first column again
D2190922_true <- tibble::rownames_to_column(D2190922_true, "ASV")

#save new dataframe of only ASV, total depth, and fish
D2190922_final <- D2190922_true[c(1,5:7)]

#190922-F2----
#turn ASV into row names, remove non-numeric columns
F2190922_true <- F2190922_1perc[,-1]
rownames(F2190922_true) <- F2190922_1perc[,1]
F2190922_true <- F2190922_true[c(2:4)]

#subset to keep only PCR reps where the total depth > 10% of the reads
F2190922_true <- F2190922_true[,colSums(F2190922_true) >= (sum(F2190922_1perc$Total.depth)*0.1)]

#get new total depth of these ASVs
F2190922_true$Total.depth <- apply(F2190922_true[,1:3],1,sum)

#add column for number of PCR reps this total depth came from
F2190922_true$Num.PCR.reps <- "3"

#add back in fish name
F2190922_true$Fish <- "190922-F2"

#make ASV the first column again
F2190922_true <- tibble::rownames_to_column(F2190922_true, "ASV")

#save new dataframe of only ASV, total depth, and fish
F2190922_final <- F2190922_true[c(1,5:7)]

#190922-G2----
#turn ASV into row names, remove non-numeric columns
G2190922_true <- G2190922_1perc[,-1]
rownames(G2190922_true) <- G2190922_1perc[,1]
G2190922_true <- G2190922_true[c(2:4)]

#subset to keep only PCR reps where the total depth > 10% of the reads
G2190922_true <- G2190922_true[,colSums(G2190922_true) >= (sum(G2190922_1perc$Total.depth)*0.1)]

#get new total depth of these ASVs
G2190922_true$Total.depth <- apply(G2190922_true[,1:3],1,sum)

#add column for number of PCR reps this total depth came from
G2190922_true$Num.PCR.reps <- "3"

#add back in fish name
G2190922_true$Fish <- "190922-G2"

#make ASV the first column again
G2190922_true <- tibble::rownames_to_column(G2190922_true, "ASV")

#save new dataframe of only ASV, total depth, and fish
G2190922_final <- G2190922_true[c(1,5:7)]

#190222-H1----
#turn ASV into row names, remove non-numeric columns
H1190222_true <- H1190222_1perc[,-1]
rownames(H1190222_true) <- H1190222_1perc[,1]
H1190222_true <- H1190222_true[c(2:4)]

#subset to keep only PCR reps where the total depth > 10% of the reads
H1190222_true <- H1190222_true[,colSums(H1190222_true) >= (sum(H1190222_1perc$Total.depth)*0.1)]

#keep only non-zero ASVs (some ASVs only present in PCR rep 3, which did not contain
#at least 1% of the total depth of the fish)
H1190222_true <- subset(H1190222_true, PCR.1>0 & PCR.2>0)

#get new total depth of these ASVs
H1190222_true$Total.depth <- apply(H1190222_true[,1:2],1,sum)

#add column for number of PCR reps this total depth came from
H1190222_true$Num.PCR.reps <- "2"

#add back in fish name
H1190222_true$Fish <- "190222-H1"

#make ASV the first column again
H1190222_true <- tibble::rownames_to_column(H1190222_true, "ASV")

#save new dataframe of only ASV, total depth, and fish
H1190222_final <- H1190222_true[c(1,4:6)]

#190922-H2----
#turn ASV into row names, remove non-numeric columns
H2190922_true <- H2190922_1perc[,-1]
rownames(H2190922_true) <- H2190922_1perc[,1]
H2190922_true <- H2190922_true[c(2:4)]

#subset to keep only PCR reps where the total depth > 10% of the reads
H2190922_true <- H2190922_true[,colSums(H2190922_true) >= (sum(H2190922_1perc$Total.depth)*0.1)]

#get new total depth of these ASVs
H2190922_true$Total.depth <- apply(H2190922_true[,1:3],1,sum)

#add column for number of PCR reps this total depth came from
H2190922_true$Num.PCR.reps <- "3"

#add back in fish name
H2190922_true$Fish <- "190922-H2"

#make ASV the first column again
H2190922_true <- tibble::rownames_to_column(H2190922_true, "ASV")

#save new dataframe of only ASV, total depth, and fish
H2190922_final <- H2190922_true[c(1,5:7)]

#190922-J2----
#turn ASV into row names, remove non-numeric columns
J2190922_true <- J2190922_1perc[,-1]
rownames(J2190922_true) <- J2190922_1perc[,1]
J2190922_true <- J2190922_true[c(2:4)]

#subset to keep only PCR reps where the total depth > 10% of the reads
J2190922_true <- J2190922_true[,colSums(J2190922_true) >= (sum(J2190922_1perc$Total.depth)*0.1)]

#get new total depth of these ASVs
J2190922_true$Total.depth <- apply(J2190922_true[,1:3],1,sum)

#add column for number of PCR reps this total depth came from
J2190922_true$Num.PCR.reps <- "3"

#add back in fish name
J2190922_true$Fish <- "190922-J2"

#make ASV the first column again
J2190922_true <- tibble::rownames_to_column(J2190922_true, "ASV")

#save new dataframe of only ASV, total depth, and fish
J2190922_final <- J2190922_true[c(1,5:7)]

#190222-M1----
#turn ASV into row names, remove non-numeric columns
M1190222_true <- M1190222_1perc[,-1]
rownames(M1190222_true) <- M1190222_1perc[,1]
M1190222_true <- M1190222_true[c(2:4)]

#subset to keep only PCR reps where the total depth > 10% of the reads
M1190222_true <- M1190222_true[,colSums(M1190222_true) >= (sum(M1190222_1perc$Total.depth)*0.1)]

#get new total depth of these ASVs
M1190222_true$Total.depth <- apply(M1190222_true[,1:3],1,sum)

#add column for number of PCR reps this total depth came from
M1190222_true$Num.PCR.reps <- "3"

#add back in fish name
M1190222_true$Fish <- "190222-M1"

#make ASV the first column again
M1190222_true <- tibble::rownames_to_column(M1190222_true, "ASV")

#save new dataframe of only ASV, total depth, and fish
M1190222_final <- M1190222_true[c(1,5:7)]

#NTC----
#turn ASV into row names, remove non-numeric columns
NTC_true <- NTC_1perc[,-1]
rownames(NTC_true) <- NTC_1perc[,1]
NTC_true <- NTC_true[c(2:4)]

#subset to keep only PCR reps where the total depth > 10% of the reads
NTC_true <- NTC_true[,colSums(NTC_true) >= (sum(NTC_1perc$Total.depth)*0.1)]

#keep only non-zero ASVs (some ASVs only present in PCR rep 1, which did not contain
#at least 1% of the total depth of the fish)
NTC_true <- subset(NTC_true, PCR.2>0 & PCR.3>0)

#get new total depth of these ASVs
NTC_true$Total.depth <- apply(NTC_true[,1:2],1,sum)

#add column for number of PCR reps this total depth came from
NTC_true$Num.PCR.reps <- "2"

#add back in fish name
NTC_true$Fish <- "NTC"

#make ASV the first column again
NTC_true <- tibble::rownames_to_column(NTC_true, "ASV")

#save new dataframe of only ASV, total depth, and fish
NTC_final <- NTC_true[c(1,4:6)]




#join all together----
all_1perc_true <- rbind(A1190222_final,A2190222_final, H1190222_final,M1190222_final,
                        A1190918_final,A3190918_final,A4190918_final,C4190918_final,
                        D4190918_final, E7190918_final,F2190918_final,B2190922_final,
                        D2190922_final, F2190922_final,G2190922_final,H2190922_final,
                        J2190922_final,NTC_final)


#add sequences
#read in data (from previously made .csv of UNIX seqs - moved into current working directory)
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
library(tidyverse)
all_1perc_true_seqs <- all_1perc_true %>% 
  inner_join(UNIXseqs, by = c("ASV" = "ASV"))

#move sequence to the front
all_1perc_true_seqs <- all_1perc_true_seqs %>% relocate(Sequence, .after = ASV)

#save output for later analysis
write_csv(all_1perc_true_seqs, "organized_ASVs/1%_true_ASVs.csv")


#reassess H1-190222, B2-190922, and A1-190918 fish----
#(to ensure retained ASVs are still present in >=1% total depth of the remaining PCR reps)
H1190222_true_1perc_again <- subset(H1190222_true, PCR.1>=(sum(H1190222_true$PCR.1)*0.01) | PCR.2>=(sum(H1190222_true$PCR.2)*0.01))

#save new dataframe of only ASV, total depth, and fish
H1190222_final_1perc_again <- H1190222_true_1perc_again[c(1,4:6)]


A1190918_true_1perc_again <- subset(A1190918_true, PCR.1>=(sum(A1190918_true$PCR.1)*0.01) | PCR.3>=(sum(A1190918_true$PCR.3)*0.01))
A1190918_final_1perc_again <- A1190918_true_1perc_again[c(1,4:6)]

B2190922_true_1perc_again <- subset(B2190922_true, PCR.1>=(sum(B2190922_true$PCR.1)*0.01) | PCR.2>=(sum(B2190922_true$PCR.2)*0.01))
B2190922_final_1perc_again <- B2190922_true_1perc_again[c(1,4:6)]



all_1perc_true_final <- rbind(A1190222_final,A2190222_final, H1190222_final_1perc_again,M1190222_final,
                              A1190918_final_1perc_again,A3190918_final,A4190918_final,C4190918_final,
                        D4190918_final, E7190918_final,F2190918_final,B2190922_final_1perc_again,
                        D2190922_final, F2190922_final,G2190922_final,H2190922_final,
                        J2190922_final,NTC_final)

#add sequences back in
#read in data (from previously made .csv of UNIX seqs - moved into current working directory)
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
library(tidyverse)
all_1perc_true_final_seqs <- all_1perc_true_final %>% 
  inner_join(UNIXseqs, by = c("ASV" = "ASV"))

#move sequence to the front
all_1perc_true_final_seqs <- all_1perc_true_final_seqs %>% relocate(Sequence, .after = ASV)

#save output for later analysis
write_csv(all_1perc_true_final_seqs, "organized_ASVs/1%_true_1%_ASVs.csv")



