#analysis of the ASVs returned from the UNIX pipeline when run with default parameters except:
#sequence length restricted to 70-105bp (as determined from multiqc report)
#phred score 26 enforced during trimming
#sequences only taken through the dereplication stage

#Previously parse through true ASVs to find only those with depths >=10% per fish
#set working directory----
getwd()

#"/Users/samanthabeal"
setwd("Documents/MSc/Bioinformatics")

ASVs <- read.csv("UNIX-opt-new/ASVdepth_1percent.csv")
fish <- ASVs[c(1,5)]
fish$ASV <- as.factor(fish$ASV)
fish$Fish <- as.factor(fish$Fish)

uniquefish <- as.data.frame(unique(fish$Fish))
names(uniquefish)[1] <- "Fish"


#subset into a dataframe for each ASV

#idea:
#something like making a data frame for each ASV (^)
#then something like if each fish is present, put yes, if not put no
#will need to have the code output to a place i can keep track of it


X <- split(ASVs, ASVs$ASV)
str(X)
Y <- lapply(seq_along(X), function(x) as.data.frame(X[[x]])[,1:5])
uniq1 <- Y[[1]]
uniq10 <- Y[[2]]
uniq11 <- Y[[3]]
uniq12 <- Y[[4]]
uniq13 <- Y[[5]]
uniq14 <- Y[[6]]
uniq15 <- Y[[7]]
uniq17 <- Y[[8]]
uniq18 <- Y[[9]]
uniq19 <- Y[[10]]
uniq2 <- Y[[11]]
uniq20 <- Y[[12]]
uniq21 <- Y[[13]]
uniq24 <- Y[[14]]
uniq25 <- Y[[15]]
uniq26 <- Y[[16]]
uniq27 <- Y[[17]]
uniq3 <- Y[[18]]
uniq4 <- Y[[19]]
uniq5 <- Y[[20]]
uniq6 <- Y[[21]]
uniq64 <- Y[[22]]
uniq7 <- Y[[23]]
uniq8 <- Y[[24]]
uniq9 <- Y[[25]]
uniq92 <- Y[[26]]



#check if fish in in each ASV df
uniq1_fish <- as.data.frame(uniquefish$Fish %in% uniq1$Fish)
names(uniq1_fish)[1] <- "uniq1"

uniq10_fish <- as.data.frame(uniquefish$Fish %in% uniq10$Fish)
names(uniq10_fish)[1] <- "uniq10"

uniq11_fish <- as.data.frame(uniquefish$Fish %in% uniq11$Fish)
names(uniq11_fish)[1] <- "uniq11"

uniq12_fish <- as.data.frame(uniquefish$Fish %in% uniq12$Fish)
names(uniq12_fish)[1] <- "uniq12"

uniq13_fish <- as.data.frame(uniquefish$Fish %in% uniq13$Fish)
names(uniq13_fish)[1] <- "uniq13"

uniq14_fish <- as.data.frame(uniquefish$Fish %in% uniq14$Fish)
names(uniq14_fish)[1] <- "uniq14"

uniq15_fish <- as.data.frame(uniquefish$Fish %in% uniq15$Fish)
names(uniq15_fish)[1] <- "uniq15"

uniq17_fish <- as.data.frame(uniquefish$Fish %in% uniq17$Fish)
names(uniq17_fish)[1] <- "uniq17"

uniq18_fish <- as.data.frame(uniquefish$Fish %in% uniq18$Fish)
names(uniq18_fish)[1] <- "uniq18"

uniq19_fish <- as.data.frame(uniquefish$Fish %in% uniq19$Fish)
names(uniq19_fish)[1] <- "uniq19"

uniq2_fish <- as.data.frame(uniquefish$Fish %in% uniq2$Fish)
names(uniq2_fish)[1] <- "uniq2"

uniq20_fish <- as.data.frame(uniquefish$Fish %in% uniq20$Fish)
names(uniq20_fish)[1] <- "uniq20"

uniq21_fish<- as.data.frame(uniquefish$Fish %in% uniq21$Fish)
names(uniq21_fish)[1] <- "uniq21"

uniq24_fish <- as.data.frame(uniquefish$Fish %in% uniq24$Fish)
names(uniq24_fish)[1] <- "uniq24"

uniq25_fish <- as.data.frame(uniquefish$Fish %in% uniq25$Fish)
names(uniq25_fish)[1] <- "uniq25"

uniq26_fish <- as.data.frame(uniquefish$Fish %in% uniq26$Fish)
names(uniq26_fish)[1] <- "uniq26"

uniq27_fish <- as.data.frame(uniquefish$Fish %in% uniq27$Fish)
names(uniq27_fish)[1] <- "uniq27"

uniq3_fish <- as.data.frame(uniquefish$Fish %in% uniq3$Fish)
names(uniq3_fish)[1] <- "uniq3"

uniq4_fish <- as.data.frame(uniquefish$Fish %in% uniq4$Fish)
names(uniq4_fish)[1] <- "uniq4"

uniq5_fish <- as.data.frame(uniquefish$Fish %in% uniq5$Fish)
names(uniq5_fish)[1] <- "uniq5"

uniq6_fish <- as.data.frame(uniquefish$Fish %in% uniq6$Fish)
names(uniq6_fish)[1] <- "uniq6"

uniq64_fish <- as.data.frame(uniquefish$Fish %in% uniq64$Fish)
names(uniq64_fish)[1] <- "uniq64"

uniq7_fish <- as.data.frame(uniquefish$Fish %in% uniq7$Fish)
names(uniq7_fish)[1] <- "uniq7"

uniq8_fish <- as.data.frame(uniquefish$Fish %in% uniq8$Fish)
names(uniq8_fish)[1] <- "uniq8"

uniq9_fish <- as.data.frame(uniquefish$Fish %in% uniq9$Fish)
names(uniq9_fish)[1] <- "uniq9"

uniq92_fish <- as.data.frame(uniquefish$Fish %in% uniq92$Fish)
names(uniq92_fish)[1] <- "uniq92"


fishperASV <- cbind(uniq1_fish,uniq10_fish,uniq11_fish,uniq12_fish,uniq13_fish,
                    uniq14_fish,uniq15_fish,uniq17_fish,uniq18_fish,uniq19_fish,
                    uniq2_fish,uniq20_fish,uniq21_fish,uniq24_fish,uniq25_fish,
                    uniq26_fish,uniq27_fish,uniq3_fish,uniq4_fish,uniq5_fish,
                    uniq6_fish,uniq64_fish,uniq7_fish,uniq8_fish,uniq9_fish,
                    uniq92_fish)


row.names(fishperASV) <- uniquefish$Fish

#change from true/false to 1/0
fishperASV <- as.numeric(fishperASV)
fishperASV <- as.integer(as.logical(fishperASV))

fishperASVmatrix <- data.matrix(fishperASV)

#plot----
library(ggplot2)
library(reshape2)

# Convert wide to long
matrix.long <- melt(fishperASVmatrix)
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
        text = element_text(size = 20),
        axis.text.x = element_text(angle=90)) +
  xlab("ASV") + ylab("Fish") +
  ggtitle("Consistently detected ASVs present in individuals at rates >= 1% total read depth") 




#compare these seqs to dada2----
#load dada2 >=1% depth seqs
dada <-  read.csv("dada2-optimization/o10_non151_ASVdepth_1percent.csv")
uniqASV <- as.data.frame(unique(ASVs[c("ASV","Sequence")]))

UNIXvdada <- as.data.frame(uniqASV$Sequence %in% dada$Sequence)
row.names(UNIXvdada) <- uniqASV$Sequence
names(UNIXvdada)[1] <- "UNIX ASV in dada?"

library(dplyr)
UNIXvdada <- tibble::rownames_to_column(UNIXvdada, "Sequence")

comp <- uniqASV %>% 
  inner_join(UNIXvdada, by = c("Sequence" = "Sequence"))



dadavUNIX<- as.data.frame(uniqASV$Sequence %in% dada$Sequence)
row.names(dadavUNIX) <- uniqASV$Sequence
names(dadavUNIX)[1] <- "dada ASV in UNIX?"

dadavUNIX <- tibble::rownames_to_column(dadavUNIX, "Sequence")
dada2 <- dada %>% 
  inner_join(dadavUNIX, by = c("Sequence" = "Sequence"))
