#ASV analysis of SINE UNIX optimization


#set wd to "/Users/samanthabeal/Documents/MSc/Bioinformatics"
getwd()
#"/Users/samanthabeal"

setwd("Documents/MSc/Bioinformatics")

#load data----
library(readr)
new <-readr::read_tsv("UNIX-opt-new/output1/ASV/derep.tsv", show_col_types = FALSE)
#rename columns to only be fish and index
names(new) = gsub(pattern = "_.*", replacement = "", x = names(new))

#pull out individual fish----
#190222-A1----
A1190222n1 <- new[c(1:4)]
names(A1190222n1)[1] <- "ASV"
names(A1190222n1)[2] <- "PCR.1"
names(A1190222n1)[3] <- "PCR.2"
names(A1190222n1)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
A1190222n1$Consistently_detected <- ifelse(A1190222n1$PCR.1==0 | A1190222n1$PCR.2== 0 | A1190222n1$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
A1190222n1_true <- subset(A1190222n1, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
A1190222n1_true$Fish <- "A1190222"
A1190222n1_true$Output <- "derep only"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
A1190222n1sum <- as.data.frame(colSums(A1190222n1[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A1190222n1sum <- as.data.frame(t(A1190222n1sum))

#add column of avg depth across all PCR reps
A1190222n1sum$Avg.depth <-apply(A1190222n1sum,1,mean)
#add columns of total depth
A1190222n1sum$Total.depth <-apply(A1190222n1sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
A1190222n1sum$Num.ASVs <- sum(A1190222n1$PCR.1>0 | A1190222n1$PCR.2>0 | A1190222n1$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
A1190222n1sum$Output <- "derep only"

#add column of if in all 3 reps
A1190222n1sum$All.reps <- "No"

#only true ASVs
A1190222n1_truesum <- as.data.frame(colSums(A1190222n1_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A1190222n1_truesum <- as.data.frame(t(A1190222n1_truesum))

#add column of avg depth across all PCR reps
A1190222n1_truesum $Avg.depth <-apply(A1190222n1_truesum,1,mean)

#add columns of total depth
A1190222n1_truesum $Total.depth <-apply(A1190222n1_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
A1190222n1_truesum$Num.ASVs <- sum(A1190222n1_true$PCR.1>0 | A1190222n1_true$PCR.2>0 | A1190222n1_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
A1190222n1_truesum$Output <- "derep only"

#add column of if in all 3 reps
A1190222n1_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
A1190222_new <- rbind(A1190222n1_truesum, A1190222n1sum)

#add column of fish ID
A1190222_new$Fish <- "190922-A1"

#190222-A2----
A2190222n1 <- new[c(1,5:7)]
names(A2190222n1)[1] <- "ASV"
names(A2190222n1)[2] <- "PCR.1"
names(A2190222n1)[3] <- "PCR.2"
names(A2190222n1)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
A2190222n1$Consistently_detected <- ifelse(A2190222n1$PCR.1==0 | A2190222n1$PCR.2== 0 | A2190222n1$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
A2190222n1_true <- subset(A2190222n1, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
A2190222n1_true$Fish <- "A2190222"
A2190222n1_true$Output <- "derep only"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
A2190222n1sum <- as.data.frame(colSums(A2190222n1[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A2190222n1sum <- as.data.frame(t(A2190222n1sum))

#add column of avg depth across all PCR reps
A2190222n1sum$Avg.depth <-apply(A2190222n1sum,1,mean)
#add columns of total depth
A2190222n1sum$Total.depth <-apply(A2190222n1sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
A2190222n1sum$Num.ASVs <- sum(A2190222n1$PCR.1>0 | A2190222n1$PCR.2>0 | A2190222n1$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
A2190222n1sum$Output <- "derep only"

#add column of if in all 3 reps
A2190222n1sum$All.reps <- "No"

#only true ASVs
A2190222n1_truesum <- as.data.frame(colSums(A2190222n1_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A2190222n1_truesum <- as.data.frame(t(A2190222n1_truesum))

#add column of avg depth across all PCR reps
A2190222n1_truesum $Avg.depth <-apply(A2190222n1_truesum,1,mean)

#add columns of total depth
A2190222n1_truesum $Total.depth <-apply(A2190222n1_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
A2190222n1_truesum$Num.ASVs <- sum(A2190222n1_true$PCR.1>0 | A2190222n1_true$PCR.2>0 | A2190222n1_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
A2190222n1_truesum$Output <- "derep only"

#add column of if in all 3 reps
A2190222n1_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
A2190222_new <- rbind(A2190222n1_truesum, A2190222n1sum)

#add column of fish ID
A2190222_new$Fish <- "190922-A2"

#190222-H1----
H1190222n1 <- new[c(1,8:10)]
names(H1190222n1)[1] <- "ASV"
names(H1190222n1)[2] <- "PCR.1"
names(H1190222n1)[3] <- "PCR.2"
names(H1190222n1)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
H1190222n1$Consistently_detected <- ifelse(H1190222n1$PCR.1==0 | H1190222n1$PCR.2== 0 | H1190222n1$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
H1190222n1_true <- subset(H1190222n1, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
H1190222n1_true$Fish <- "H1190222"
H1190222n1_true$Output <- "derep only"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
H1190222n1sum <- as.data.frame(colSums(H1190222n1[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
H1190222n1sum <- as.data.frame(t(H1190222n1sum))

#add column of avg depth across all PCR reps
H1190222n1sum$Avg.depth <-apply(H1190222n1sum,1,mean)
#add columns of total depth
H1190222n1sum$Total.depth <-apply(H1190222n1sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
H1190222n1sum$Num.ASVs <- sum(H1190222n1$PCR.1>0 | H1190222n1$PCR.2>0 | H1190222n1$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
H1190222n1sum$Output <- "derep only"

#add column of if in all 3 reps
H1190222n1sum$All.reps <- "No"

#only true ASVs
H1190222n1_truesum <- as.data.frame(colSums(H1190222n1_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
H1190222n1_truesum <- as.data.frame(t(H1190222n1_truesum))

#add column of avg depth across all PCR reps
H1190222n1_truesum $Avg.depth <-apply(H1190222n1_truesum,1,mean)

#add columns of total depth
H1190222n1_truesum $Total.depth <-apply(H1190222n1_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
H1190222n1_truesum$Num.ASVs <- sum(H1190222n1_true$PCR.1>0 | H1190222n1_true$PCR.2>0 | H1190222n1_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
H1190222n1_truesum$Output <- "derep only"

#add column of if in all 3 reps
H1190222n1_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
H1190222_new <- rbind(H1190222n1_truesum, H1190222n1sum)

#add column of fish ID
H1190222_new$Fish <- "190922-H1"

#190222-M1----
M1190222n1 <- new[c(1,11:13)]
names(M1190222n1)[1] <- "ASV"
names(M1190222n1)[2] <- "PCR.1"
names(M1190222n1)[3] <- "PCR.2"
names(M1190222n1)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
M1190222n1$Consistently_detected <- ifelse(M1190222n1$PCR.1==0 | M1190222n1$PCR.2== 0 | M1190222n1$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
M1190222n1_true <- subset(M1190222n1, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
M1190222n1_true$Fish <- "M1190222"
M1190222n1_true$Output <- "derep only"


#compare size of two dfs (# ASVs and total reads)
#all ASVs
M1190222n1sum <- as.data.frame(colSums(M1190222n1[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
M1190222n1sum <- as.data.frame(t(M1190222n1sum))

#add column of avg depth across all PCR reps
M1190222n1sum$Avg.depth <-apply(M1190222n1sum,1,mean)
#add columns of total depth
M1190222n1sum$Total.depth <-apply(M1190222n1sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
M1190222n1sum$Num.ASVs <- sum(M1190222n1$PCR.1>0 | M1190222n1$PCR.2>0 | M1190222n1$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
M1190222n1sum$Output <- "derep only"

#add column of if in all 3 reps
M1190222n1sum$All.reps <- "No"

#only true ASVs
M1190222n1_truesum <- as.data.frame(colSums(M1190222n1_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
M1190222n1_truesum <- as.data.frame(t(M1190222n1_truesum))

#add column of avg depth across all PCR reps
M1190222n1_truesum $Avg.depth <-apply(M1190222n1_truesum,1,mean)

#add columns of total depth
M1190222n1_truesum $Total.depth <-apply(M1190222n1_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
M1190222n1_truesum$Num.ASVs <- sum(M1190222n1_true$PCR.1>0 | M1190222n1_true$PCR.2>0 | M1190222n1_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
M1190222n1_truesum$Output <- "derep only"

#add column of if in all 3 reps
M1190222n1_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
M1190222_new <- rbind(M1190222n1_truesum, M1190222n1sum)

#add column of fish ID
M1190222_new$Fish <- "190922-M1"

#190918-A1----
A1190912n1 <- new[c(1,14:16)]
names(A1190912n1)[1] <- "ASV"
names(A1190912n1)[2] <- "PCR.1"
names(A1190912n1)[3] <- "PCR.2"
names(A1190912n1)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
A1190912n1$Consistently_detected <- ifelse(A1190912n1$PCR.1==0 | A1190912n1$PCR.2== 0 | A1190912n1$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
A1190912n1_true <- subset(A1190912n1, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
A1190912n1_true$Fish <- "A1190918"
A1190912n1_true$Output <- "derep only"
#compare size of two dfs (# ASVs and total reads)
#all ASVs
A1190912n1sum <- as.data.frame(colSums(A1190912n1[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A1190912n1sum <- as.data.frame(t(A1190912n1sum))

#add column of avg depth across all PCR reps
A1190912n1sum$Avg.depth <-apply(A1190912n1sum,1,mean)
#add columns of total depth
A1190912n1sum$Total.depth <-apply(A1190912n1sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
A1190912n1sum$Num.ASVs <- sum(A1190912n1$PCR.1>0 | A1190912n1$PCR.2>0 | A1190912n1$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
A1190912n1sum$Output <- "derep only"

#add column of if in all 3 reps
A1190912n1sum$All.reps <- "No"

#only true ASVs
A1190912n1_truesum <- as.data.frame(colSums(A1190912n1_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A1190912n1_truesum <- as.data.frame(t(A1190912n1_truesum))

#add column of avg depth across all PCR reps
A1190912n1_truesum $Avg.depth <-apply(A1190912n1_truesum,1,mean)

#add columns of total depth
A1190912n1_truesum $Total.depth <-apply(A1190912n1_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
A1190912n1_truesum$Num.ASVs <- sum(A1190912n1_true$PCR.1>0 | A1190912n1_true$PCR.2>0 | A1190912n1_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
A1190912n1_truesum$Output <- "derep only"

#add column of if in all 3 reps
A1190912n1_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
A1190918_new <- rbind(A1190912n1_truesum, A1190912n1sum)

#add column of fish ID
A1190918_new$Fish <- "190918-A1"

#190918-A3----
A3190912n1 <- new[c(1,17:19)]
names(A3190912n1)[1] <- "ASV"
names(A3190912n1)[2] <- "PCR.1"
names(A3190912n1)[3] <- "PCR.2"
names(A3190912n1)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
A3190912n1$Consistently_detected <- ifelse(A3190912n1$PCR.1==0 | A3190912n1$PCR.2== 0 | A3190912n1$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
A3190912n1_true <- subset(A3190912n1, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
A3190912n1_true$Fish <- "A3190918"
A3190912n1_true$Output <- "derep only"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
A3190912n1sum <- as.data.frame(colSums(A3190912n1[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A3190912n1sum <- as.data.frame(t(A3190912n1sum))

#add column of avg depth across all PCR reps
A3190912n1sum$Avg.depth <-apply(A3190912n1sum,1,mean)
#add columns of total depth
A3190912n1sum$Total.depth <-apply(A3190912n1sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
A3190912n1sum$Num.ASVs <- sum(A3190912n1$PCR.1>0 | A3190912n1$PCR.2>0 | A3190912n1$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
A3190912n1sum$Output <- "derep only"

#add column of if in all 3 reps
A3190912n1sum$All.reps <- "No"

#only true ASVs
A3190912n1_truesum <- as.data.frame(colSums(A3190912n1_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A3190912n1_truesum <- as.data.frame(t(A3190912n1_truesum))

#add column of avg depth across all PCR reps
A3190912n1_truesum $Avg.depth <-apply(A3190912n1_truesum,1,mean)

#add columns of total depth
A3190912n1_truesum $Total.depth <-apply(A3190912n1_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
A3190912n1_truesum$Num.ASVs <- sum(A3190912n1_true$PCR.1>0 | A3190912n1_true$PCR.2>0 | A3190912n1_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
A3190912n1_truesum$Output <- "derep only"

#add column of if in all 3 reps
A3190912n1_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
A3190918_new <- rbind(A3190912n1_truesum, A3190912n1sum)

#add column of fish ID
A3190918_new$Fish <- "190918-A3"

#190918-A4----
A4190912n1 <- new[c(1,20:22)]
names(A4190912n1)[1] <- "ASV"
names(A4190912n1)[2] <- "PCR.1"
names(A4190912n1)[3] <- "PCR.2"
names(A4190912n1)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
A4190912n1$Consistently_detected <- ifelse(A4190912n1$PCR.1==0 | A4190912n1$PCR.2== 0 | A4190912n1$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
A4190912n1_true <- subset(A4190912n1, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
A4190912n1_true$Fish <- "A4190918"
A4190912n1_true$Output <- "derep only"


#compare size of two dfs (# ASVs and total reads)
#all ASVs
A4190912n1sum <- as.data.frame(colSums(A4190912n1[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A4190912n1sum <- as.data.frame(t(A4190912n1sum))

#add column of avg depth across all PCR reps
A4190912n1sum$Avg.depth <-apply(A4190912n1sum,1,mean)
#add columns of total depth
A4190912n1sum$Total.depth <-apply(A4190912n1sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
A4190912n1sum$Num.ASVs <- sum(A4190912n1$PCR.1>0 | A4190912n1$PCR.2>0 | A4190912n1$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
A4190912n1sum$Output <- "derep only"

#add column of if in all 3 reps
A4190912n1sum$All.reps <- "No"

#only true ASVs
A4190912n1_truesum <- as.data.frame(colSums(A4190912n1_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A4190912n1_truesum <- as.data.frame(t(A4190912n1_truesum))

#add column of avg depth across all PCR reps
A4190912n1_truesum $Avg.depth <-apply(A4190912n1_truesum,1,mean)

#add columns of total depth
A4190912n1_truesum $Total.depth <-apply(A4190912n1_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
A4190912n1_truesum$Num.ASVs <- sum(A4190912n1_true$PCR.1>0 | A4190912n1_true$PCR.2>0 | A4190912n1_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
A4190912n1_truesum$Output <- "derep only"

#add column of if in all 3 reps
A4190912n1_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
A4190918_new <- rbind(A4190912n1_truesum, A4190912n1sum)

#add column of fish ID
A4190918_new$Fish <- "190918-A4"

#190918-C4----
C4190912n1 <- new[c(1,23:25)]
names(C4190912n1)[1] <- "ASV"
names(C4190912n1)[2] <- "PCR.1"
names(C4190912n1)[3] <- "PCR.2"
names(C4190912n1)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
C4190912n1$Consistently_detected <- ifelse(C4190912n1$PCR.1==0 | C4190912n1$PCR.2== 0 | C4190912n1$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
C4190912n1_true <- subset(C4190912n1, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
C4190912n1_true$Fish <- "C4190918"
C4190912n1_true$Output <- "derep only"


#compare size of two dfs (# ASVs and total reads)
#all ASVs
C4190912n1sum <- as.data.frame(colSums(C4190912n1[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
C4190912n1sum <- as.data.frame(t(C4190912n1sum))

#add column of avg depth across all PCR reps
C4190912n1sum$Avg.depth <-apply(C4190912n1sum,1,mean)
#add columns of total depth
C4190912n1sum$Total.depth <-apply(C4190912n1sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
C4190912n1sum$Num.ASVs <- sum(C4190912n1$PCR.1>0 | C4190912n1$PCR.2>0 | C4190912n1$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
C4190912n1sum$Output <- "derep only"

#add column of if in all 3 reps
C4190912n1sum$All.reps <- "No"

#only true ASVs
C4190912n1_truesum <- as.data.frame(colSums(C4190912n1_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
C4190912n1_truesum <- as.data.frame(t(C4190912n1_truesum))

#add column of avg depth across all PCR reps
C4190912n1_truesum $Avg.depth <-apply(C4190912n1_truesum,1,mean)

#add columns of total depth
C4190912n1_truesum $Total.depth <-apply(C4190912n1_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
C4190912n1_truesum$Num.ASVs <- sum(C4190912n1_true$PCR.1>0 | C4190912n1_true$PCR.2>0 | C4190912n1_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
C4190912n1_truesum$Output <- "derep only"

#add column of if in all 3 reps
C4190912n1_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
C4190918_new <- rbind(C4190912n1_truesum, C4190912n1sum)

#add column of fish ID
C4190918_new$Fish <- "190918-C4"

#190918-D4----
D4190912n1 <- new[c(1,26:28)]
names(D4190912n1)[1] <- "ASV"
names(D4190912n1)[2] <- "PCR.1"
names(D4190912n1)[3] <- "PCR.2"
names(D4190912n1)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
D4190912n1$Consistently_detected <- ifelse(D4190912n1$PCR.1==0 | D4190912n1$PCR.2== 0 | D4190912n1$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
D4190912n1_true <- subset(D4190912n1, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
D4190912n1_true$Fish <- "D4190918"
D4190912n1_true$Output <- "derep only"


#compare size of two dfs (# ASVs and total reads)
#all ASVs
D4190912n1sum <- as.data.frame(colSums(D4190912n1[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
D4190912n1sum <- as.data.frame(t(D4190912n1sum))

#add column of avg depth across all PCR reps
D4190912n1sum$Avg.depth <-apply(D4190912n1sum,1,mean)
#add columns of total depth
D4190912n1sum$Total.depth <-apply(D4190912n1sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
D4190912n1sum$Num.ASVs <- sum(D4190912n1$PCR.1>0 | D4190912n1$PCR.2>0 | D4190912n1$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
D4190912n1sum$Output <- "derep only"

#add column of if in all 3 reps
D4190912n1sum$All.reps <- "No"

#only true ASVs
D4190912n1_truesum <- as.data.frame(colSums(D4190912n1_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
D4190912n1_truesum <- as.data.frame(t(D4190912n1_truesum))

#add column of avg depth across all PCR reps
D4190912n1_truesum$Avg.depth <-apply(D4190912n1_truesum,1,mean)

#add columns of total depth
D4190912n1_truesum$Total.depth <-apply(D4190912n1_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
D4190912n1_truesum$Num.ASVs <- sum(D4190912n1_true$PCR.1>0 | D4190912n1_true$PCR.2>0 | D4190912n1_true$PCR.3>0)

#add column of output (so can plot with UNIX for comparison)
D4190912n1_truesum$Output <- "derep only"

#add column of if in all 3 reps
D4190912n1_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
D4190918_new <- rbind(D4190912n1_truesum, D4190912n1sum)

#add column of fish ID
D4190918_new$Fish <- "190918-D4"

#190918-E7----
E7190912n1 <- new[c(1,29:31)]
names(E7190912n1)[1] <- "ASV"
names(E7190912n1)[2] <- "PCR.1"
names(E7190912n1)[3] <- "PCR.2"
names(E7190912n1)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
E7190912n1$Consistently_detected <- ifelse(E7190912n1$PCR.1==0 | E7190912n1$PCR.2== 0 | E7190912n1$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
E7190912n1_true <- subset(E7190912n1, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
E7190912n1_true$Fish <- "E7190918"
E7190912n1_true$Output <- "derep only"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
E7190912n1sum <- as.data.frame(colSums(E7190912n1[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
E7190912n1sum <- as.data.frame(t(E7190912n1sum))

#add column of avg depth across all PCR reps
E7190912n1sum$Avg.depth <-apply(E7190912n1sum,1,mean)
#add columns of total depth
E7190912n1sum$Total.depth <-apply(E7190912n1sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
E7190912n1sum$Num.ASVs <- sum(E7190912n1$PCR.1>0 | E7190912n1$PCR.2>0 | E7190912n1$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
E7190912n1sum$Output <- "derep only"

#add column of if in all 3 reps
E7190912n1sum$All.reps <- "No"

#only true ASVs
E7190912n1_truesum <- as.data.frame(colSums(E7190912n1_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
E7190912n1_truesum <- as.data.frame(t(E7190912n1_truesum))

#add column of avg depth across all PCR reps
E7190912n1_truesum $Avg.depth <-apply(E7190912n1_truesum,1,mean)

#add columns of total depth
E7190912n1_truesum $Total.depth <-apply(E7190912n1_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
E7190912n1_truesum$Num.ASVs <- sum(E7190912n1_true$PCR.1>0 | E7190912n1_true$PCR.2>0 | E7190912n1_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
E7190912n1_truesum$Output <- "derep only"

#add column of if in all 3 reps
E7190912n1_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
E7190918_new <- rbind(E7190912n1_truesum, E7190912n1sum)

#add column of fish ID
E7190918_new$Fish <- "190918-E7"

#190918-F2----
F2190912n1 <- new[c(1,32:34)]
names(F2190912n1)[1] <- "ASV"
names(F2190912n1)[2] <- "PCR.1"
names(F2190912n1)[3] <- "PCR.2"
names(F2190912n1)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
F2190912n1$Consistently_detected <- ifelse(F2190912n1$PCR.1==0 | F2190912n1$PCR.2== 0 | F2190912n1$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
F2190912n1_true <- subset(F2190912n1, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
F2190912n1_true$Fish <- "F2190918"
F2190912n1_true$Output <- "derep only"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
F2190912n1sum <- as.data.frame(colSums(F2190912n1[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
F2190912n1sum <- as.data.frame(t(F2190912n1sum))

#add column of avg depth across all PCR reps
F2190912n1sum$Avg.depth <-apply(F2190912n1sum,1,mean)
#add columns of total depth
F2190912n1sum$Total.depth <-apply(F2190912n1sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
F2190912n1sum$Num.ASVs <- sum(F2190912n1$PCR.1>0 | F2190912n1$PCR.2>0 | F2190912n1$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
F2190912n1sum$Output <- "derep only"

#add column of if in all 3 reps
F2190912n1sum$All.reps <- "No"

#only true ASVs
F2190912n1_truesum <- as.data.frame(colSums(F2190912n1_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
F2190912n1_truesum <- as.data.frame(t(F2190912n1_truesum))

#add column of avg depth across all PCR reps
F2190912n1_truesum $Avg.depth <-apply(F2190912n1_truesum,1,mean)

#add columns of total depth
F2190912n1_truesum $Total.depth <-apply(F2190912n1_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
F2190912n1_truesum$Num.ASVs <- sum(F2190912n1_true$PCR.1>0 | F2190912n1_true$PCR.2>0 | F2190912n1_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
F2190912n1_truesum$Output <- "derep only"

#add column of if in all 3 reps
F2190912n1_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
F2190918_new <- rbind(F2190912n1_truesum, F2190912n1sum)

#add column of fish ID
F2190918_new$Fish <- "190918-F2"

#190922-B2----
B2190922n1 <- new[c(1,35:37)]
names(B2190922n1)[1] <- "ASV"
names(B2190922n1)[2] <- "PCR.1"
names(B2190922n1)[3] <- "PCR.2"
names(B2190922n1)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
B2190922n1$Consistently_detected <- ifelse(B2190922n1$PCR.1==0 | B2190922n1$PCR.2== 0 | B2190922n1$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
B2190922n1_true <- subset(B2190922n1, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
B2190922n1_true$Fish <- "B2190922"
B2190922n1_true$Output <- "derep only"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
B2190922n1sum <- as.data.frame(colSums(B2190922n1[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
B2190922n1sum <- as.data.frame(t(B2190922n1sum))

#add column of avg depth across all PCR reps
B2190922n1sum$Avg.depth <-apply(B2190922n1sum,1,mean)
#add columns of total depth
B2190922n1sum$Total.depth <-apply(B2190922n1sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
B2190922n1sum$Num.ASVs <- sum(B2190922n1$PCR.1>0 | B2190922n1$PCR.2>0 | B2190922n1$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
B2190922n1sum$Output <- "derep only"

#add column of if in all 3 reps
B2190922n1sum$All.reps <- "No"

#only true ASVs
B2190922n1_truesum <- as.data.frame(colSums(B2190922n1_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
B2190922n1_truesum <- as.data.frame(t(B2190922n1_truesum))

#add column of avg depth across all PCR reps
B2190922n1_truesum $Avg.depth <-apply(B2190922n1_truesum,1,mean)

#add columns of total depth
B2190922n1_truesum $Total.depth <-apply(B2190922n1_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
B2190922n1_truesum$Num.ASVs <- sum(B2190922n1_true$PCR.1>0 | B2190922n1_true$PCR.2>0 | B2190922n1_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
B2190922n1_truesum$Output <- "derep only"

#add column of if in all 3 reps
B2190922n1_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
B2190922_new <- rbind(B2190922n1_truesum, B2190922n1sum)

#add column of fish ID
B2190922_new$Fish <- "190922-B2"

#190922-D2----
D2190922n1 <- new[c(1,38:40)]
names(D2190922n1)[1] <- "ASV"
names(D2190922n1)[2] <- "PCR.1"
names(D2190922n1)[3] <- "PCR.2"
names(D2190922n1)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
D2190922n1$Consistently_detected <- ifelse(D2190922n1$PCR.1==0 | D2190922n1$PCR.2== 0 | D2190922n1$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
D2190922n1_true <- subset(D2190922n1, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
D2190922n1_true$Fish <- "D2190922"
D2190922n1_true$Output <- "derep only"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
D2190922n1sum <- as.data.frame(colSums(D2190922n1[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
D2190922n1sum <- as.data.frame(t(D2190922n1sum))

#add column of avg depth across all PCR reps
D2190922n1sum$Avg.depth <-apply(D2190922n1sum,1,mean)
#add columns of total depth
D2190922n1sum$Total.depth <-apply(D2190922n1sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
D2190922n1sum$Num.ASVs <- sum(D2190922n1$PCR.1>0 | D2190922n1$PCR.2>0 | D2190922n1$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
D2190922n1sum$Output <- "derep only"

#add column of if in all 3 reps
D2190922n1sum$All.reps <- "No"

#only true ASVs
D2190922n1_truesum <- as.data.frame(colSums(D2190922n1_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
D2190922n1_truesum <- as.data.frame(t(D2190922n1_truesum))

#add column of avg depth across all PCR reps
D2190922n1_truesum $Avg.depth <-apply(D2190922n1_truesum,1,mean)

#add columns of total depth
D2190922n1_truesum $Total.depth <-apply(D2190922n1_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
D2190922n1_truesum$Num.ASVs <- sum(D2190922n1_true$PCR.1>0 | D2190922n1_true$PCR.2>0 | D2190922n1_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
D2190922n1_truesum$Output <- "derep only"

#add column of if in all 3 reps
D2190922n1_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
D2190922_new <- rbind(D2190922n1_truesum, D2190922n1sum)

#add column of fish ID
D2190922_new$Fish <- "190922-D2"

#190922-F2----
F2190922n1 <- new[c(1,41:43)]
names(F2190922n1)[1] <- "ASV"
names(F2190922n1)[2] <- "PCR.1"
names(F2190922n1)[3] <- "PCR.2"
names(F2190922n1)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
F2190922n1$Consistently_detected <- ifelse(F2190922n1$PCR.1==0 | F2190922n1$PCR.2== 0 | F2190922n1$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
F2190922n1_true <- subset(F2190922n1, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
F2190922n1_true$Fish <- "F2190922"
F2190922n1_true$Output <- "derep only"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
F2190922n1sum <- as.data.frame(colSums(F2190922n1[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
F2190922n1sum <- as.data.frame(t(F2190922n1sum))

#add column of avg depth across all PCR reps
F2190922n1sum$Avg.depth <-apply(F2190922n1sum,1,mean)
#add columns of total depth
F2190922n1sum$Total.depth <-apply(F2190922n1sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
F2190922n1sum$Num.ASVs <- sum(F2190922n1$PCR.1>0 | F2190922n1$PCR.2>0 | F2190922n1$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
F2190922n1sum$Output <- "derep only"

#add column of if in all 3 reps
F2190922n1sum$All.reps <- "No"

#only true ASVs
F2190922n1_truesum <- as.data.frame(colSums(F2190922n1_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
F2190922n1_truesum <- as.data.frame(t(F2190922n1_truesum))

#add column of avg depth across all PCR reps
F2190922n1_truesum $Avg.depth <-apply(F2190922n1_truesum,1,mean)

#add columns of total depth
F2190922n1_truesum $Total.depth <-apply(F2190922n1_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
F2190922n1_truesum$Num.ASVs <- sum(F2190922n1_true$PCR.1>0 | F2190922n1_true$PCR.2>0 | F2190922n1_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
F2190922n1_truesum$Output <- "derep only"

#add column of if in all 3 reps
F2190922n1_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
F2190922_new <- rbind(F2190922n1_truesum, F2190922n1sum)

#add column of fish ID
F2190922_new$Fish <- "190922-F2"

#190922-G2----
G2190922n1 <- new[c(1,44:46)]
names(G2190922n1)[1] <- "ASV"
names(G2190922n1)[2] <- "PCR.1"
names(G2190922n1)[3] <- "PCR.2"
names(G2190922n1)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
G2190922n1$Consistently_detected <- ifelse(G2190922n1$PCR.1==0 | G2190922n1$PCR.2== 0 | G2190922n1$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
G2190922n1_true <- subset(G2190922n1, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
G2190922n1_true$Fish <- "G2190922"
G2190922n1_true$Output <- "derep only"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
G2190922n1sum <- as.data.frame(colSums(G2190922n1[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
G2190922n1sum <- as.data.frame(t(G2190922n1sum))

#add column of avg depth across all PCR reps
G2190922n1sum$Avg.depth <-apply(G2190922n1sum,1,mean)
#add columns of total depth
G2190922n1sum$Total.depth <-apply(G2190922n1sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
G2190922n1sum$Num.ASVs <- sum(G2190922n1$PCR.1>0 | G2190922n1$PCR.2>0 | G2190922n1$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
G2190922n1sum$Output <- "derep only"

#add column of if in all 3 reps
G2190922n1sum$All.reps <- "No"

#only true ASVs
G2190922n1_truesum <- as.data.frame(colSums(G2190922n1_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
G2190922n1_truesum <- as.data.frame(t(G2190922n1_truesum))

#add column of avg depth across all PCR reps
G2190922n1_truesum $Avg.depth <-apply(G2190922n1_truesum,1,mean)

#add columns of total depth
G2190922n1_truesum $Total.depth <-apply(G2190922n1_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
G2190922n1_truesum$Num.ASVs <- sum(G2190922n1_true$PCR.1>0 | G2190922n1_true$PCR.2>0 | G2190922n1_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
G2190922n1_truesum$Output <- "derep only"

#add column of if in all 3 reps
G2190922n1_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
G2190922_new <- rbind(G2190922n1_truesum, G2190922n1sum)

#add column of fish ID
G2190922_new$Fish <- "190922-G2"

#190922-H2----
H2190922n1 <- new[c(1,47:49)]
names(H2190922n1)[1] <- "ASV"
names(H2190922n1)[2] <- "PCR.1"
names(H2190922n1)[3] <- "PCR.2"
names(H2190922n1)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
H2190922n1$Consistently_detected <- ifelse(H2190922n1$PCR.1==0 | H2190922n1$PCR.2== 0 | H2190922n1$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
H2190922n1_true <- subset(H2190922n1, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
H2190922n1_true$Fish <- "H2190922"
H2190922n1_true$Output <- "derep only"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
H2190922n1sum <- as.data.frame(colSums(H2190922n1[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
H2190922n1sum <- as.data.frame(t(H2190922n1sum))

#add column of avg depth across all PCR reps
H2190922n1sum$Avg.depth <-apply(H2190922n1sum,1,mean)
#add columns of total depth
H2190922n1sum$Total.depth <-apply(H2190922n1sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
H2190922n1sum$Num.ASVs <- sum(H2190922n1$PCR.1>0 | H2190922n1$PCR.2>0 | H2190922n1$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
H2190922n1sum$Output <- "derep only"

#add column of if in all 3 reps
H2190922n1sum$All.reps <- "No"

#only true ASVs
H2190922n1_truesum <- as.data.frame(colSums(H2190922n1_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
H2190922n1_truesum <- as.data.frame(t(H2190922n1_truesum))

#add column of avg depth across all PCR reps
H2190922n1_truesum $Avg.depth <-apply(H2190922n1_truesum,1,mean)

#add columns of total depth
H2190922n1_truesum $Total.depth <-apply(H2190922n1_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
H2190922n1_truesum$Num.ASVs <- sum(H2190922n1_true$PCR.1>0 | H2190922n1_true$PCR.2>0 | H2190922n1_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
H2190922n1_truesum$Output <- "derep only"

#add column of if in all 3 reps
H2190922n1_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
H2190922_new <- rbind(H2190922n1_truesum, H2190922n1sum)

#add column of fish ID
H2190922_new$Fish <- "190922-H2"

#190922-J2----
J2190922n1 <- new[c(1,50:52)]
names(J2190922n1)[1] <- "ASV"
names(J2190922n1)[2] <- "PCR.1"
names(J2190922n1)[3] <- "PCR.2"
names(J2190922n1)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
J2190922n1$Consistently_detected <- ifelse(J2190922n1$PCR.1==0 | J2190922n1$PCR.2== 0 | J2190922n1$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
J2190922n1_true <- subset(J2190922n1, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
J2190922n1_true$Fish <- "J2190922"
J2190922n1_true$Output <- "derep only"


#compare size of two dfs (# ASVs and total reads)
#all ASVs
J2190922n1sum <- as.data.frame(colSums(J2190922n1[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
J2190922n1sum <- as.data.frame(t(J2190922n1sum))

#add column of avg depth across all PCR reps
J2190922n1sum$Avg.depth <-apply(J2190922n1sum,1,mean)
#add columns of total depth
J2190922n1sum$Total.depth <-apply(J2190922n1sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
J2190922n1sum$Num.ASVs <- sum(J2190922n1$PCR.1>0 | J2190922n1$PCR.2>0 | J2190922n1$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
J2190922n1sum$Output <- "derep only"

#add column of if in all 3 reps
J2190922n1sum$All.reps <- "No"

#only true ASVs
J2190922n1_truesum <- as.data.frame(colSums(J2190922n1_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
J2190922n1_truesum <- as.data.frame(t(J2190922n1_truesum))

#add column of avg depth across all PCR reps
J2190922n1_truesum $Avg.depth <-apply(J2190922n1_truesum,1,mean)

#add columns of total depth
J2190922n1_truesum $Total.depth <-apply(J2190922n1_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
J2190922n1_truesum$Num.ASVs <- sum(J2190922n1_true$PCR.1>0 | J2190922n1_true$PCR.2>0 | J2190922n1_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
J2190922n1_truesum$Output <- "derep only"

#add column of if in all 3 reps
J2190922n1_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
J2190922_new <- rbind(J2190922n1_truesum, J2190922n1sum)

#add column of fish ID
J2190922_new$Fish <- "190922-J2"






#NTC----
NTCn1 <- new[c(1,53:55)]
names(NTCn1)[1] <- "ASV"
names(NTCn1)[2] <- "PCR.1"
names(NTCn1)[3] <- "PCR.2"
names(NTCn1)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
NTCn1$Consistently_detected <- ifelse(NTCn1$PCR.1==0 | NTCn1$PCR.2== 0 | NTCn1$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
NTCn1_true <- subset(NTCn1, Consistently_detected == "true",
                     select=c("ASV","PCR.1","PCR.2","PCR.3"))
NTCn1_true$Fish <- "NTC"
NTCn1_true$Output <- "derep only"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
NTCn1sum <- as.data.frame(colSums(NTCn1[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
NTCn1sum <- as.data.frame(t(NTCn1sum))

#add column of avg depth across all PCR reps
NTCn1sum$Avg.depth <-apply(NTCn1sum,1,mean)

#add columns of total depth
NTCn1sum$Total.depth <-apply(NTCn1sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
NTCn1sum$Num.ASVs <- sum(NTCn1$PCR.1>0 | NTCn1$PCR.2>0 | NTCn1$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
NTCn1sum$Output <- "new"

#add column of if in all 3 reps
NTCn1sum$All.reps <- "No"

#only true ASVs
NTCn1_truesum <- as.data.frame(colSums(NTCn1_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
NTCn1_truesum <- as.data.frame(t(NTCn1_truesum))

#add column of avg depth across all PCR reps
NTCn1_truesum $Avg.depth <-apply(NTCn1_truesum,1,mean)

#add columns of total depth
NTCn1_truesum $Total.depth <-apply(NTCn1_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
NTCn1_truesum$Num.ASVs <- sum(NTCn1_true$PCR.1>0 | NTCn1_true$PCR.2>0 | NTCn1_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
NTCn1_truesum$Output <- "new"

#add column of if in all 3 reps
NTCn1_truesum$All.reps <- "Yes"

#NTC add data frames together
NTC_new <- rbind(NTCn1_truesum, NTCn1sum)

#add column of fish ID
NTC_new$Fish <- "NTC"

#put it all all together ----
all_new <- rbind(A1190222_new,A2190222_new,H1190222_new,M1190222_new,
                   A1190918_new,A3190918_new,A4190918_new,C4190918_new,
                   D4190918_new,E7190918_new,F2190918_new,B2190922_new,
                   D2190922_new, F2190922_new,G2190922_new,H2190922_new,
                 J2190922_new,NTC_new)

#reorganizae
rownames(all_new)<-c(1:36)
library(tidyverse)
all_new <- all_new %>% relocate(Fish, .before = PCR.1)

#save for future analysis
write_csv(all_new, "UNIX-opt-new/ASVanalysis.csv")


#change data types
all_new$Fish <- as.factor(all_new$Fish)
all_new$Output <- as.factor(all_new$Output)
all_new$All.reps <- as.factor(all_new$All.reps)


#plot with bars for individual fish
#want to add a label to each bar showing how many ASVs this average was taken from
avgASV <- aggregate(Num.ASVs ~ Fish + All.reps, data = all_new, FUN = mean)

avgdepth <- aggregate(Avg.depth ~ Fish + All.reps, data = all_new, FUN = mean)


#plot with bars for individual fish
ggplot(avgdepth, aes(fill=All.reps, y=Avg.depth, x=Fish)) + 
  geom_bar(position="dodge", stat="identity") + 
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Pipeline") + ylab("Average read depth") + 
  labs(fill= "ASV detected in all PCR replicates?") +
  geom_text(aes(label = round(avgASV$Num.ASVs, digits = 0)), position=position_dodge(width=0.9), vjust=-0.25, size=3)

#average across ALL fish
avgallASV <- aggregate(Num.ASVs ~ All.reps, data = all_new, FUN = mean)

avgalldepth <- aggregate(Avg.depth ~ All.reps, data = all_new, FUN = mean)

ggplot(avgalldepth, aes(fill=All.reps, y=Avg.depth, x=All.reps)) + 
  geom_bar(position="dodge", stat="identity") + 
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("ASV detected in all PCR replicates?") + ylab("Average read depth") + 
  labs(fill= "ASV detected in all PCR replicates?") +
  ggtitle("ASV analysis of new") +
  guides(fill="none") +
  geom_text(aes(label = round(avgallASV$Num.ASVs, digits = 0)), position=position_dodge(width=0.9), vjust=-0.25, size=3)



#only looking at consistently detected ASVs----
alltruen1 <- rbind(A1190222n1_true,A2190222n1_true,H1190222n1_true,M1190222n1_true,
                   A1190912n1_true,A3190912n1_true,A4190912n1_true,C4190912n1_true,
                   D4190912n1_true,E7190912n1_true,F2190912n1_true,B2190922n1_true,
                   D2190922n1_true,F2190922n1_true,G2190922n1_true,H2190922n1_true,
                   J2190922n1_true,NTCn1_true)

#save for future analysis
write_csv(alltruen1, "UNIX-opt-new/TrueASVs.csv")


uniqueASVt1 <- unique(alltruen1$ASV)
length(unique(alltruen1$ASV))
#8891


#add column of total depth across all PCR reps
alltruen1$Total.depth <-apply(alltruen1[,2:4],1,sum)
sum(alltruen1$Total.depth)

#total depth of unique ASVs per technical replicate summed across all individuals
#569,344

#if want to see total depth per fish:
aggregate(Total.depth ~ Fish, data = alltruen1, FUN = sum)







