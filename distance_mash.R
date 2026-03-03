setwd(dir="~/marina/julius-pt/Julius_jan/run3/pt3_julius/all")
distances <- read.table(file="~/marina/julius-pt/Julius_jan/run3/pt3_julius/all/mash.txt")
distances <- distances[!distances$V3>0.05,]
#distances[!duplicated(distances$V2),]
library(tidyr)
distances <- separate(data=distances, col ="V5", into = c("shared_mashes", "sketch_size"))
distances <- distances %>%
  arrange(V2, desc(shared_mashes))

#up first arrange according to sample and shared mashes then do duplicates 
top_hits <- distances[!duplicated(distances$V2),] # top hits per genome 
write.table(top_hits,file="tophits.txt")

##file from efetch

spc <- read.csv(file="~/Downloads/speciator (15).csv")
#names red in libre office <- "Refname", "Accession", "Species"

ref <- read.csv(file="~/marina/julius-pt/Julius_jan/run3/pt3_julius/all/ref.txt", sep="\t")
top_hits2 <- read.table(file="tophits.txt",header = TRUE)
colnames(top_hits2) <- c("Ref","Sample","Mash_Distance","p_value","shared_mashes","sketch_size")

library(dplyr)
colnames(ref) 
colnames(top_hits2)
top_hits2 <- left_join(top_hits2,ref,join_by("Ref"=="Refname")) ##name the mash hits from entrez

sts <- read.table(file="~/marina/julius-pt/Julius_jan/run3/pt3_julius/all/mlst.txt",header = TRUE)
top_hits2$Sample <- gsub(pattern = ".fasta",replacement = "",x=top_hits2$Sample)

top_hits3 <- left_join(top_hits2,spc[,c(2,4)],join_by("Sample"=="Genome.Name"))
colnames(sts)

sts$Sample <- gsub(pattern = ".fasta",replacement = "",sts$Sample)

top_hits4 <- left_join(top_hits3,sts,by="Sample")

write.table(top_hits4,file="tophits4_final.txt")



all_amr <- read.csv(file="~/marina/julius-pt/Julius_jan/run3/pt3_julius/all/amr_finder/all.txt",
                    sep = "\t", header=TRUE)

all_amr <- all_amr[!all_amr$Name=="Name",]

all_amr <- all_amr[!all_amr$Name=="Name",]

write.csv(all_amr,file="amr_genes.csv", row.names = FALSE)

all_amr$Subclass %>% unique()







