
setwd(dir ="~/marina/prof_dam/new_results_kleb")

sam_names <- read.csv(file= "~/marina/prof_dam/kleb/kleb_list2.txt", sep = "\t", header=TRUE)

sts <- read.csv(file="~/marina/prof_dam/kleb/kleborate_results/klebsiella_pneumo_complex_output.txt",
                sep = "\t", header = TRUE)


kleborate <- sts
species <- kleborate$species
spc_freq <- summary(as.factor(species))
spc_freq <- as.data.frame(spc_freq)
colnames(spc_freq) <- "freq"
rnames <- row.names(spc_freq)
spc_freq2 <- as.data.frame(cbind(rnames, spc_freq$freq))
kala <- c("#00CED1", "#FC8D62", "#8DA0CB", "#E78AC3")
rnames[c(1,2,3,4)] <- c("K.pneumoniae", "K.q.similipneumoniae", "K.v.tropica", "K.v.variicola")
#par(mar = c(5, 4, 4, 2))  # Bottom, Left, Top, Right
pie(x=spc_freq$freq,labels = spc_freq$freq, main="Pie Chart showing species distribution", col=kala)
legend("right",rnames, cex = 0.5,fill =kala, 
       bty = "n",       # Remove box border (optional)
       inset = c(0.5, 0.1),    # Slightly inset the legend from the margin
       x.intersp = 0.3, # Reduce horizontal spacing between boxes and text
       y.intersp = 0.5) #Reduce vertical spacing between legend items





library(dplyr)



all_amr <- read.csv(file="~/marina/prof_dam/kleb/amrfinder_two/all.txt",
                    sep = "\t", header=TRUE)

all_amr <- all_amr[!all_amr$Name=="Name",]

all_amr <- all_amr[!all_amr$Name=="Name",]
all_amr$Subclass %>% unique()
all_amr$Subclass[all_amr$Subclass=="SPECTINOMYCIN/STREPTOMYCIN"] <- "SPE/STR"
all_amr$Subclass[all_amr$Subclass=="AMIKACIN/KANAMYCIN/QUINOLONE/TOBRAMYCIN"] <- "AMI/KAN/QUI/TOB"
all_amr$Subclass[all_amr$Subclass=="CLINDAMYCIN/ERYTHROMYCIN/STREPTOGRAMIN B"] <- "CLI/ERY/STRG_B"
all_amr$Subclass[all_amr$Subclass=="AZITHROMYCIN/ERYTHROMYCIN/SPIRAMYCIN/TELITHROMYCIN"] <- "AZI/ERY/SPI/TEL"
all_amr$Subclass[all_amr$Subclass=="AZITHROMYCIN/ERYTHROMYCIN/STREPTOGRAMIN"] <- "AZI/ERY/STRG"
all_amr$Subclass[all_amr$Subclass=="PHENICOL/QUINOLONE"] <- "PHE/QUI"
all_amr$Subclass[all_amr$Subclass=="AZTREONAM/CEFIDEROCOL/CEPHALOSPORIN"] <- "AZI/CEFD/CEP"
all_amr$Subclass[all_amr$Subclass=="AMPICILLIN/CHLORAMPHENICOL/QUINOLONE/RIFAMPIN/TETRACYCLINE"] <- "AMP/CHL/QUI/RIF/TET"
all_amr$Subclass[all_amr$Subclass=="GENTAMICIN/KANAMYCIN/TOBRAMYCIN"] <- "GEN/KAN/TOB"
all_amr$Subclass[all_amr$Subclass=="AMIKACIN/GENTAMICIN/KANAMYCIN/TOBRAMYCIN"] <- "AMI/GEN/KAN/TOB"
#all_amr$Subclass[all_amr$Subclass=="PHENICOL/QUINOLONE"] <- "PHE/QUI

all_amr$Class %>% unique()

all_amr$Class[all_amr$Class=="LINCOSAMIDE/MACROLIDE/STREPTOGRAMIN"] <- "LIN/MAC/STRG"
all_amr$Class[all_amr$Class=="AMINOGLYCOSIDE/QUINOLONE"] <- "AMG/QUI"
all_amr$Class[all_amr$Class=="PHENICOL/QUINOLONE"] <- "PHE/QUI"
all_amr$Class[all_amr$Class=="MACROLIDE/STREPTOGRAMIN"] <- "MAC/STRG"






all_point <- all_amr[all_amr$Subtype=="POINT",]
genes_amr <- all_amr[all_amr$Subtype=="AMR",]
genes_amr$Class %>% table()

#number of amr genes

unique(genes_amr$Element.symbol) %>% length()
unique(genes_amr$Class) %>% length()

noG_sample <- genes_amr$Name %>% table() %>% as.data.frame()
noG_sample <- noG_sample[order(noG_sample$Freq, decreasing = TRUE),]
summary(noG_sample$Freq)

amr_genes_freq <- genes_amr$Element.symbol %>% table() %>% as.data.frame()
amr_genes_freq <- amr_genes_freq[order(amr_genes_freq$Freq, decreasing = TRUE),]

write.csv(file="amr_genes_freq.csv",amr_genes_freq, row.names = FALSE)

genes <- unique(genes_amr$Element.symbol)
samples <- unique(genes_amr$Name)


##presence or absence matrix for for genes
p_b_matrix <- matrix(nrow = length(samples), ncol = length(genes)+1) %>% as.data.frame()
#p_b_matrix <- matrix(nrow = length(samples), ncol = length(genes)) %>% as.data.frame()
p_b_matrix$V1 <- samples
colnames(p_b_matrix)[1] <- "sample"
colnames(p_b_matrix)[-1] <- genes
#colnames(p_b_matrix) <- genes
for (sample in samples){
  df <- genes_amr[genes_amr$Name==sample,]
  ind <- colnames(p_b_matrix) %in% df$Element.symbol %>% which()
  p_b_matrix[p_b_matrix$sample==sample, ind] <- 1
  ind <- append(ind,values =1,after = length(ind)) #prevents sample name from being over written with 0
  #print(ind)
  p_b_matrix[p_b_matrix$sample==sample, -ind] <- 0
  
  
}




genes_p_b <- as.data.frame(p_b_matrix)

##checks for accuracy

k <- apply(genes_p_b[,-1],2,as.numeric)
apply(k,2,sum) #frequency should be the same as amr_genes_freq

#genes_p_b <- as.numeric(genes_p_b)

#genes_p_b[,-1] <- lapply(genes_p_b[,-1], function(x) as.numeric(as.character(x)))
gn <- apply(genes_p_b[,-1],2,as.numeric)
genes_p_b[,-1] <-  gn


ord <- apply(genes_p_b[,-1],2,sum) %>% as.data.frame()
ind_inc <- rownames(ord)[order(ord$.,decreasing = TRUE)]

ind_inc <- append(ind_inc,value="sample",after=0) ## 0 dds be4 first element

gene_p_b <- genes_p_b[,ind_inc] 

#sample arrangement
#rownames(gene_p_b) <- gene_p_b$sample
ord_sample <- apply(genes_p_b[,-1],1,sum) %>% as.data.frame()
ind_inc_sample <- rownames(ord_sample)[order(ord_sample$.,decreasing = TRUE)]
ind_inc_sample <- as.numeric(ind_inc_sample)
ind_inc_sample <- gene_p_b$sample[ind_inc_sample]
library(ggplot2)
library(reshape)

#melt(gene_p_b)

df_long <- melt(gene_p_b,id.vars = "sample", variable.name = "Gene", value.name = "Presence")
#df_long$Gene <- factor(df_long$Gene,levels =ind_inc )
colnames(df_long) <- c("sample","Gene","Presence")
df_long$sample <- factor(df_long$sample,levels =ind_inc_sample )

k <- genes_amr[,c("Class","Element.symbol","Subclass")]

antiB_class <- k[!duplicated(k$Element.symbol),] 

df_long$antibiotic <- NA
df_long$subclass <- NA
for (i in 1:nrow(df_long)){
  ind <- (antiB_class$Element.symbol %in% df_long$Gene[i]) %>% which
  df_long$antibiotic[i] <- antiB_class$Class[ind]
  df_long$subclass[i] <- antiB_class$Subclass[ind]
}







df_present <- df_long %>%
  filter(Presence == 1)

##custom color palete

pub_palette_20 <- c(
  # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#CC79A7", # reddish purple
  "#999999", # grey
  "#117733", # dark green
  # teal
  "#88CCEE", # light blue
  "#DDCC77", # sand
  "#AA4499", # purple-magenta
  # wine red
  "#6699CC", # steel blue
  "#661100", # brown
  "#888800", # olive
  "#FFAABB", # light pink
  "#66CCEE", # cyan
  "#222222"  # near black
)

# Stacked barplot by antibiotic class
ggplot(df_present, aes(x = sample, fill = antibiotic)) +
  geom_bar() +
  labs(x = "Sample",
       y = "Number of Resistance Genes",
       fill = "Antibiotic Class") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


pub_palette_18 <- c(
  # Original Set3
  "#8DD3C7", "#377EB8", "#BEBADA", "#FB8072", "#80B1D3",
  "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
  "#CCEBC5", "#1B9E77",
  # 6 extra distinct colors
  "#E41A1C", # strong red
  "#1B9E77", # dark teal
  "#FF7F00", # vivid orange
  "#984EA3", # purple
  "#377EB8", # steel blue
  "#A65628"  # brown
)



stacked_bar <- ggplot(df_present[df_present$Presence == 1, ], 
       aes(x = sample, fill = antibiotic)) +
  geom_bar(position = "stack") +
  scale_fill_manual(values = pub_palette_18) +  # Publication-friendly palette
  theme_minimal(base_size = 20) +
  labs(
    x = "Sample",
    y = "Number of Resistance Genes",
    fill = "Antibiotic Class",
    title = "Antimicrobial Resistance Genes Across Samples"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "bold"),
    panel.grid.major.x = element_blank()
  )

#ggsave("amr_bar_gene.png", plot = stacked_bar, width = 10, height = 6, dpi = 600)
ggsave(
  filename = "amr_bar_gene.svg",
  plot =stacked_bar,
  device = "svg",
  width = 14,
  height = 10
  
)


library(cowplot)
stacked_bar_nolgnd <- stacked_bar + theme(legend.position = "none")
ggsave(
  filename = "stacked_bar_nolgnd.svg",
  plot = stacked_bar_nolgnd,
  device = "svg",
  width = 14,
  height = 10
  
)

stacked_barlegend_only <- get_legend(stacked_bar)
ggsave("stacked_bar.svg",stacked_barlegend_only , width = 10, height = 5)


#sort the df according to antibiotic class
df_long <- arrange(df_long, antibiotic)
df_long$Gene <- factor(df_long$Gene, levels = unique(df_long$Gene))


#head(df_long)
htmp <- ggplot(df_long, aes(x = sample, y = Gene)) +
  geom_tile(aes(alpha = Presence,fill=antibiotic),color="white") +
  scale_alpha(range = c(0, 1), guide = "none") + # fade out absent values
  scale_x_discrete(expand = c(0, 0)) + 
  theme_minimal(base_size = 25) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.5, size = 12, face= "bold"),
    axis.title = element_blank(),
    legend.title = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 12, face= "bold")
  ) 
#ggsave("htmap_gene_profiles.png", plot = htmp, width = 10, height = 6, dpi = 600)

library(cowplot)
htmp_nolgnd <- htmp + theme(legend.position = "none")
ggsave(
  filename = "htmp_nolgnd.svg",
  plot = htmp_nolgnd,
  device = "svg",
  width = 14,
  height = 12
  
)

htmp_horizontal <- htmp +
  theme(
    legend.position = "bottom",   # place legend below the plot
    legend.direction = "horizontal",  # make it horizontal
    legend.text =  element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 16)
  )



legend_onlyhtmp <- get_legend(htmp_horizontal)
ggsave("htmplgnd.svg", legend_onlyhtmp, width = 20, height = 5)






##presence or absence matrix for antibiotics
anti_B <- unique(genes_amr$Class)
samples <- unique(genes_amr$Name)

p_b_matrix <- matrix(nrow = length(samples), ncol = length(anti_B)+1) %>% as.data.frame()
#p_b_matrix <- matrix(nrow = length(samples), ncol = length(genes)) %>% as.data.frame()
p_b_matrix$V1 <- samples
colnames(p_b_matrix)[1] <- "sample"
colnames(p_b_matrix)[-1] <- anti_B
#colnames(p_b_matrix) <- genes
for (sample in samples){
  df <- genes_amr[genes_amr$Name==sample,]
  ind <- colnames(p_b_matrix) %in% df$Class %>% which()
  p_b_matrix[p_b_matrix$sample==sample, ind] <- 1
  ind <- append(ind,values =1,after = length(ind)) #prevents sample name from being over written with 0
  #print(ind)
  p_b_matrix[p_b_matrix$sample==sample, -ind] <- 0
  
  
}

anti_B_p_b <- p_b_matrix
colnames(anti_B_p_b ) %>% as.data.frame()

m <- apply(anti_B_p_b[,-1],2,as.numeric)

apply(m,2,sum) %>% as.data.frame() ##frequencies of the diffent antibiotic classes

ord_antiBclass <- apply(anti_B_p_b[-1],2,sum)  %>% as.data.frame() 

ord_antiBclass <- rownames(ord_antiBclass)[order(ord_antiBclass$.,decreasing = TRUE)]
anti_B_p_b <-  cbind(anti_B_p_b$sample, anti_B_p_b[,ord_antiBclass]) ## arranged according to resistance
colnames(anti_B_p_b)[1] <- "sample"
write.csv(anti_B_p_b,file="anti_B_presence.csv", row.names = FALSE)



gn <- apply(anti_B_p_b[,-1],2,as.numeric)
anti_B_p_b[,-1] <-  gn



df_long_present <- df_long[df_long$Presence==1,]
subclas <- df_long_present$subclass %>% table() %>% as.data.frame()
ord_subclas <- subclas$.[order(subclas$Freq,decreasing = TRUE)]

### creating a subset of samples (with just one etrnace per drrg)
##found samples with the multiple genes confering resistancce to the same drug, skews the burden
## want one coun per drug per sample


df_long_small <- data.frame() #only one count per drug per sample

number_uniq_res <- vector()
for (sample in unique(df_long_present$sample)){
  df_sub  <- df_long_present[df_long_present$sample==sample,]
  df_sub <- df_sub[!duplicated(df_sub$subclass),]
  number_uniq_res <- append(number_uniq_res,nrow(df_sub), length(number_uniq_res))
  df_long_small <- rbind(df_long_small,df_sub)
}






##geom tile for family





df_long$subclass <- factor(df_long$subclass, levels = ord_subclas)
 
ggplot(df_long, aes(x = sample, y = subclass)) +
  geom_tile(aes(alpha = Presence),fill = "red", color="white") +
  scale_alpha(range = c(0, 1), guide = "none") + # fade out absent values
  scale_x_discrete(expand = c(0, 0)) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.5, size = 6, face= "bold"),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 8, face= "bold")
  ) 




df_long_present <- df_long_small[df_long_small$Presence==1,]
subclas <- df_long_present$subclass %>% table() %>% as.data.frame()
ord_subclas <- subclas$.[order(subclas$Freq,decreasing = TRUE)]

df_long$subclass <- factor(df_long$subclass, levels = ord_subclas)

plt_subclass <- ggplot(df_long, aes(x = sample, y = subclass)) +
  geom_tile(aes(alpha = Presence),fill = "red", color="white") +
  scale_alpha(range = c(0, 1), guide = "none") + # fade out absent values
  scale_x_discrete(expand = c(0, 0)) + 
  theme_minimal(base_size = 25) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.5, size = 12, face= "bold"),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 10, face= "bold")
  )+
  labs(title = "Resistance Presence or absence per individual Antibiotic")

#ggsave("burden_per_antibiotic.png", plot = plt_subclass, width = 10, height = 6, dpi = 600)
ggsave(
  filename = "burden_per_antibiotic.svg",
  plot = plt_subclass,
  device = "svg",
  width = 14,
  height = 10
  
)


drugclas_freq <- genes_amr$Class %>% table()
drugclas_freq <- as.data.frame(drugclas_freq) 
drugclas_freq$Freq %>% sort(decreasing = TRUE)
drugclas_freq <- drugclas_freq[order(drugclas_freq$Freq,decreasing = TRUE), ]
class_freq <- drugclas_freq$.

df_long2 <- left_join(df_long,drugclas_freq,join_by(antibiotic==.), relationship ="many-to-one")

df_long2 <- df_long2[df_long2$Presence==1,]
Gene_freq <- df_long2$Gene %>% table() %>% as.data.frame()


df_long_uniqueG <- df_long2[!duplicated(df_long2$Gene),]
uniqueG <- df_long_uniqueG$antibiotic %>% table() %>% as.data.frame()
colnames(uniqueG)[2] <- "Unique_Gene_Counts"

df_long2 <- left_join(df_long2,uniqueG,join_by(antibiotic==.), relationship ="many-to-one")

df_class <- left_join(df_class,uniqueG,join_by(antibiotic==.), relationship ="many-to-one")
gene_order <- df_class$antibiotic[order(df_class$Presence, decreasing = TRUE)]
df_class$antibiotic <- factor(df_class$antibiotic,levels=gene_order)

##plot for burden of resistance genes per antibiotic class

Hist_antBclass <- ggplot(df_class, aes(x = antibiotic, y = Presence, fill = antibiotic)) +
  geom_col(width = 0.7,show.legend = FALSE) +
  geom_text(
    aes(label = Unique_Gene_Counts,),
    hjust = -0.2,
    size = 3,
    fontface = "bold"
  ) +
  coord_flip() +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  ) +
  labs(
    y = "Number of Genes",
    x = "Antibiotic Class",
    title = "Burden of Resistance genes per Antibiotic Class"
  )

ggsave(
  filename = "Hist_antBclass.svg",
  plot = Hist_antBclass,
  device = "svg",
  width = 10,
  height = 10
  
)




#making  color strip for the phylogroups

phy_num <- unique(sts_sub$Phylogroup)

cols20 <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173",
  "#5254a3", "#9c9ede", "#6b6ecf", "#e7969c", "#31a354"
)

phylogroup_colors <- c(
  "Klebsiella pneumoniae"  = "#1b9e77",  # green
  "Klebsiella quasipneumoniae subsp. similipneumoniae" = "#d95f02",  # orange
  "Klebsiella variicola subsp. variicola" = "#7570b3"  # purple
)



##color strip for species  ##adopted from ecoli phyogroups too lazy to change

spc <- c("Klebsiella pneumoniae","Klebsiella quasipneumoniae subsp. similipneumoniae", "Klebsiella variicola subsp. variicola")

as.list(spc)

spc_itol <- kleborate[,c(1,2)]

spc_itol$color_phylo <- NA

for (i in 1:nrow(spc_itol)){
  ind <- (spc==spc_itol$species[i]) %>% which()
  #print(ind)
  grp <- spc[ind]
  spc_itol$color_phylo[i] <- paste(spc_itol$strain[i],phylogroup_colors[[grp]])
}

write.csv(spc_itol[,"color_phylo"], file = "color_strip_species.csv",row.names = FALSE)


### Virulence score

kleborate <- read.csv(file="~/marina/prof_dam/kleb/kleborate_results/klebsiella_pneumo_complex_output.txt",
                sep = "\t", header = TRUE)




#col_pts <- list(my_palette[1:length(points_ant)]) 

#8015 label #0000ff

y_only <- paste("label_background ","#B2FF66",sep=" " )
y_c_c <- paste("label_background ","#FF66B2",sep=" " )
aero <- paste("label_background ","#FF7F00",sep=" " )
a_y <- paste("label_background ","#0000ff",sep=" " )
all_3 <- paste("label_background ","#ff0000", sep=" " )

y_only <- paste("label_background ","rgba(128,0,128,0.5)",sep=" " )
y_c_c <- paste("label_background ","rgba(0,0,255,0.5)",sep=" " )
aero <- paste("label_background ","rgba(0,255,0,0.5)",sep=" " )
a_y <- paste("label_background ","rgba(255,165,0,0.5)",sep=" " )
all_3 <- paste("label_background ","rgba(255,0,0,0.5)", sep=" " )




sym <- c(y_only,y_c_c,aero,a_y,all_3)


names(sym) <- c(1,2,3,4,5)

#col_pts <- list(","#00ff00","#0000ff")


#names(col_pts) <- points_ant


kleborate_sub <- kleborate[kleborate$virulence_score>0,]
kleborate_sub$vr_itol <- NA


for (i in 1:nrow(kleborate_sub)){
  ind <- (as.numeric(names(sym))==kleborate_sub$virulence_score[i]) %>% which()
  itol <- sym[[ind]]
  itol_point <- paste(kleborate_sub$strain[i],itol)
  kleborate_sub$vr_itol[i] <- itol_point
}

kleborate_sub$vr_itol

write.table(kleborate_sub$vr_itol,file="virulence,score.txt",
            row.names = FALSE)

##making a binary dataset for itol 

#kleborate[kleborate$virulence_score
##its ybt, clbt, aero in that order
kleborate$v_itol <- NA
for (i in 1:dim(kleborate)[1]){
  if (kleborate$virulence_score[i]==0){
    kleborate$v_itol[i] <- paste(kleborate$strain[i],-1,-1,-1, sep = ",")
  }else if((kleborate$virulence_score[i]==1)){
    kleborate$v_itol[i] <- paste(kleborate$strain[i],1,-1,-1, sep = ",")
  }else if ((kleborate$virulence_score[i]==2)){
    kleborate$v_itol[i] <- paste(kleborate$strain[i],1,1,-1, sep = ",")
    
  }else if ((kleborate$virulence_score[i]==3)){
    kleborate$v_itol[i] <- paste(kleborate$strain[i],-1,-1,1, sep = ",")
  }else if ((kleborate$virulence_score[i]==4)){
    kleborate$v_itol[i] <- paste(kleborate$strain[i],1,-1,1, sep = ",")
  }else if ((kleborate$virulence_score[i]==5)){
    kleborate$v_itol[i] <- paste(kleborate$strain[i],1,1,1, sep = ",")
  }
}
  
write.table(kleborate$v_itol,file="virulence_for_binary.txt",
            row.names = FALSE)        
          




 ###PLASMIDS

plasmid1 <- read.csv(file="~/marina/prof_dam/kleb/plasmids_matrix.txt",
                     sep="\t", header=TRUE)
colnames(plasmid1)[1]
plasmid1$X.FILE <- gsub(".fasta","",plasmid1$X.FILE)
colnames(plasmid1)[1] <- "sample"
bar_chart <- paste(plasmid1$sample,plasmid1$NUM_FOUND,plasmid1$NUM_FOUND, sep = ",")

write.table(bar_chart,file="barc_chart_plasmids",
            row.names = FALSE)



pre_abs <- function(column){
  column[column != "."] <- 1
  column[column == "."] <- 0
  return(as.numeric(column)) # forces the return value to be numeric
  
}

V1 <- plasmid1[,c(1,2)]

plasmid1[] <- sapply(plasmid1, pre_abs)
plasmid1[,c(1,2)] <- V1



colnames(plasmid1)[-c(1,2)]



#plas_fam <- read.csv2(file="~/Downloads/plasmid_family_classification (1).csv", sep = ",", header = TRUE)





family <- c(
    "Col-like",   # Col.BS512._1
    "Col-like",   # Col.MG828._1
    "Col-like",   # Col.Ye4449._1
    "Col156",     # Col156_1
    "Col440II",   # Col440II_1
    "Col440I",    # Col440I_1
    "Col8282",    # Col8282_1
    "ColKP3",     # ColKP3_1
    "ColRNAI",    # ColRNAI_1
    "IncB/O/K/Z", # IncB.O.K.Z_1
    "IncB/O/K/Z", # IncB.O.K.Z_3
    "IncFIA",     # IncFIA.HI1._1_HI1
    "IncFIA",     # IncFIA_1
    "IncFIB",     # IncFIB.AP001918._1
    "IncFIB",     # IncFIB.K._1_Kpn3
    "IncFIB",     # IncFIB.Mar._1_pNDM.Mar
    "IncFIB",     # IncFIB.pQil._1_pQil
    "IncFIC/FII", # IncFIC.FII._1
    "IncFII",     # IncFII.pCRY._1_pCRY
    "IncFII",     # IncFII.pCTU2._1_pCTU2
    "IncFII",     # IncFII_1_pKP91
    "IncHI1B",    # IncHI1B_1_pNDM.MAR
    "IncL/M",     # IncL.M.pMU407._1_pMU407
    "IncN",       # IncN_1
    "IncR",       # IncR_1
    "IncX3",      # IncX3_1
    "IncY",       # IncY_1
    "p0111",      # p0111_1
    "pESA2",      # pESA2_1
    "repA_KPC"    # repA_1_pKPC.2
  )







# family <- c(
#   rep("Col", 9),
#   rep("IncB", 2),
#   rep("IncF", 10),
#   "IncI", rep("IncX", 3), "IncY"
# )

colnames(plasmid1)[-c(1,2)]

plasmid_df <- data.frame(
  plasmid = plasmid1,
  family = family
)

#plas_fam <- plasmid1



plas_fam <- apply(plasmid1[-c(1,2)],2,sum) %>% as.data.frame()
plas_fam$family <- family
colnames(plas_fam)[1] <- "freq"


rownames(plas_fam)

#colnames(plas_fam)

pfrefamily<- cbind(rownames(plas_fam),plas_fam$freq,plas_fam$family) %>% as.data.frame()

colnames(pfrefamily) <- c("Plasmid","Freq","family")

K <- order(as.numeric(pfrefamily$Freq))
pfrefamily[K,]
pfrefamily[] <- pfrefamily[K,]

plas_order <- pfrefamily$Plasmid

pfrefamily$Plasmid <- as.factor(pfrefamily$Plasmid)

pfrefamily$Plasmid <- factor(pfrefamily$Plasmid, levels =plas_order) #ggplot plots according to the factors and reorders
#them according to alphabetical order, so to get the ryt order, coarce
#colum of interest into a factor and specify the ordering of the levels using the order of interest.

#pfrefamily$Plasmid <- factor(pfrefamily$Plasmid, levels = pfrefamily$Plasmid)
pfrefamily$Freq <- as.numeric(pfrefamily$Freq)

write.csv2(pfrefamily,file="pfrefamily.txt")

library(ggplot2)

#plas_order <- pfrefamily$Plasmid

plasmid_plt <- ggplot(pfrefamily, aes(x = Plasmid, y = Freq, fill = family)) +
  geom_bar(stat = "identity") +
  labs(title = "Individual Plasmid Frequencies by Colored by Family",
       x = "Plasmid",
       y = "Freq",
       fill = "Plasmid Family") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(colour="red", size=10, 
                                face="bold"))


ggsave("Plasmid.svg", 
       plot = plasmid_plt, width = 10, height = 7,
       device = "svg",
       )

sts$escherichia__ectyper__Serotype %>% table()






sts_sub$ST %>% table()
sts_sub$ST %>% unique()

###point mutations

my_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#984EA3",
                "#00CED1", "#FF00FF", "#FFFF33", "#B2FF66", "#FF66B2") ##
all_point$Class %>% unique()
points_ant <- all_point$Class %>% unique()

#col_pts <- list(my_palette[1:length(points_ant)]) 

all_point$pt_muts <- NA

sym <- c(1,2,3) ##adjust depending on point antibiotics unique
sym <- list(1,2,3)
names(sym) <- points_ant

col_pts <- list("#ff0000","#00ff00","#0000ff")


names(col_pts) <- points_ant




for (i in 1:nrow(all_point)){
  ind <- (points_ant==all_point$Class[i]) %>% which()
  drug <- points_ant[ind]
  dr_sym <- sym[[drug]]
  itol_point <- paste(all_point$Name[i],dr_sym,10,col_pts[[drug]],1,0.5, sep = ",")
  all_point$pt_muts[i] <- itol_point
}


write.table(all_point$pt_muts,file="point_muts.txt",
            row.names = FALSE)
###color gradient for number of amrgenes

amr_gene_freq <- apply(gene_p_b[,-1],1,sum)







itol_color_grad <- paste(gene_p_b$sample,amr_gene_freq) %>% as.data.frame()


write.csv(file="colorstrip_number_resistance.csv",itol_color_grad, row.names = FALSE)

itol_color_grad2 <- cbind(gene_p_b$sample,amr_gene_freq) %>% as.data.frame()
itol_color_grad2$amr_range <- NA
itol_color_grad2$amr_gene_freq <- as.numeric(itol_color_grad2$amr_gene_freq)

for (i in 1:nrow(itol_color_grad2)){
  if(itol_color_grad2$amr_gene_freq[i] <= 4){
    itol_color_grad2$amr_range[i] <- 4
  } else if((itol_color_grad2$amr_gene_freq[i]) >4 & (itol_color_grad2$amr_gene_freq[i] < 10)){
    itol_color_grad2$amr_range[i] <- 9
  } else{itol_color_grad2$amr_range[i] <- 14}
}

itol_amr_nu_cluster <- paste(itol_color_grad2$V1,itol_color_grad2$amr_range) %>% as.data.frame()

write.csv(file="colorstrip_range_resistance.csv",itol_amr_nu_cluster, row.names = FALSE)


###stacked bar plot for amr point mutations

points_mut <- unique(all_point$Element.symbol)
samples <- unique(all_point$Name)


##presence or absence matrix for for genes
p_b_matrix <- matrix(nrow = length(samples), ncol = length(points_mut)+1) %>% as.data.frame()
#p_b_matrix <- matrix(nrow = length(samples), ncol = length(genes)) %>% as.data.frame()
p_b_matrix$V1 <- samples
colnames(p_b_matrix)[1] <- "sample"
colnames(p_b_matrix)[-1] <- points_mut
#colnames(p_b_matrix) <- genes
for (sample in samples){
  df <- all_point[all_point$Name==sample,]
  ind <- colnames(p_b_matrix) %in% df$Element.symbol %>% which()
  p_b_matrix[p_b_matrix$sample==sample, ind] <- 1
  ind <- append(ind,values =1,after = length(ind)) #prevents sample name from being over written with 0
  #print(ind)
  p_b_matrix[p_b_matrix$sample==sample, -ind] <- 0
  
  
}

point_p_b <- as.data.frame(p_b_matrix)


df_long <- melt(point_p_b,id.vars = "sample", variable.name = "Mutation", value.name = "Presence")
#df_long$Gene <- factor(df_long$Gene,levels =ind_inc )
colnames(df_long) <- c("sample","Mutation","Presence")

pt_fre <- apply(point_p_b[,-1],1,sum)

paste(point_p_b$sample,pt_fre)

ind_inc_sample <- point_p_b$sample[order(pt_fre,decreasing = TRUE)]

df_long$sample <- factor(df_long$sample,levels =ind_inc_sample )

k <- all_point[,c("Class","Element.symbol","Subclass")]





pt_class <- k[!duplicated(k$Element.symbol),] 

df_long$antibiotic <- NA
df_long$subclass <- NA
for (i in 1:nrow(df_long)){
  ind <- (pt_class$Element.symbol %in% df_long$Mutation[i]) %>% which
  df_long$antibiotic[i] <- pt_class$Class[ind]
  df_long$subclass[i] <- pt_class$Subclass[ind]
}

colors <- c("#654321", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9")


bar_points <- ggplot(df_long[df_long$Presence == 1, ], 
                     aes(x = sample, fill = subclass)) +
  geom_bar(position = "stack") +
  scale_fill_manual(values = colors) +  # Publication-friendly palette
  theme_minimal(base_size = 12) +
  labs(
    x = "Sample",
    y = "Number of Point Mutations",
    fill = "Antibiotic Class",
    title = "Antimicrobial Resistance Point mutations Across Samples"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,face = 'bold'),
    panel.grid.major.x = element_blank()
  )

ggsave(
  filename = "stacked_bar_points.svg",
  plot = bar_points,
  device = "svg",
  width = 10,
  height = 6
  
)



#ggsave("stacked_bar_points.png", plot = bar_points, width = 10, height = 6, dpi = 600)





##number of point mutations
pts <- all_point$Element.symbol %>% table() %>% as.data.frame()
colnames(pts)[1] <-"Element.symbol" 

k <- all_point[!duplicated(all_point$Element.symbol),c(1,7,13)]

k <- merge(pts,k,"Element.symbol")

write.table(k[,c(1,2,4)],file="pts_freq.txt")
##carbapenem

carba <- genes_amr[genes_amr$Subclass=="CARBAPENEM",]
write.table(carba[,c(1,7)] %>% as.data.frame(), file="cbp.txt", row.names = FALSE)
carba$Name %>% unique()

