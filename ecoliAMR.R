
setwd(dir ="~/marina/prof_dam/new_results_R")
species <- read.csv(file="~/marina/prof_dam/12dec_profdam/all/total_species.csv",header=TRUE)
colnames(species)
spcs <- table(species$Species.Name) %>% as.data.frame()
write.table(file="species_freq.txt", spcs, row.names = FALSE)





library(dplyr)

sam_names <- read.csv(file= "~/marina/prof_dam/ecoli/ecoli_list.txt", sep = "\t", header=TRUE)

sts <- read.csv(file="~/marina/prof_dam/ecoli/escherichia_output.txt",
                sep = "\t", header = TRUE)

##sts frequency

sts_freq <- sts$escherichia__mlst_achtman__ST %>% table() %>% as.data.frame()
sts_freq <- sts_freq[order(sts_freq$Freq,decreasing = TRUE),]
sum(sts_freq$Freq)

all_amr <- read.csv(file="~/marina/prof_dam/ecoli/amrfinder_two/all.txt",
                    sep = "\t", header=TRUE)

all_amr <- all_amr[!all_amr$Name=="Name",]

all_point <- all_amr[all_amr$Subtype=="POINT",]
genes_amr <- all_amr[all_amr$Subtype=="AMR",]
amr_genes_freq <- genes_amr$Element.symbol %>% table() %>% as.data.frame()

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
amr_freq <- apply(k,2,sum) #frequency should be the same as amr_genes_freq

write.table(file="amr_genes_freq.txt", as.data.frame(amr_freq))
drugclas_freq <- genes_amr$Class %>% table()

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
# 
# pub_palette_14 <- c(
#   "#E6194B", # red
#   "#3CB44B", # green
#   "#FFE119", # yellow
#   "#0082C8", # blue
#   "#F58231", # orange
#   "#911EB4", # purple
#   "#46F0F0", # cyan
#   "#F032E6", # magenta
#   "#D2F53C", # lime
#   "#FABEBE", # pink
#   "#008080", # teal
#   "#E6BEFF", # lavender
#   "#AA6E28", # brown
#   "#800000"  # maroon
# )
# 
# pub_palette_20 <- c(
#   "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
#   "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
#   "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173",
#   "#3182bd", "#31a354", "#756bb1", "#636363", "#e6550d"
# )
# 
# pub_14 <- c(
#   "#E41A1C", # red
#   "#377EB8", # blue
#   "#4DAF4A", # green
#   "#984EA3", # purple
#   "#FF7F00", # orange
#   "#FFFF33", # yellow
#   "#A65628", # brown
#   "#F781BF", # pink
#   "#999999", # grey
#   "#66C2A5", # teal
#   "#FFD92F", # bright yellow
#   "#E7298A", # magenta
#   "#1B9E77", # dark green
#   "#D95F02"  # dark orange
# )

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
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
  "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
  "#CCEBC5", "#FFED6F",
  # 6 extra distinct colors
  "#E41A1C", # strong red
  "#1B9E77", # dark teal
  "#FF7F00", # vivid orange
  "#984EA3", # purple
  "#377EB8", # steel blue
  "#A65628"  # brown
)



stacked_bar <- ggplot(df_long[df_long$Presence == 1, ], 
       aes(x = sample, fill = antibiotic)) +
  geom_bar(position = "stack") +
  scale_fill_manual(values = pub_palette_18) +  # Publication-friendly palette
  theme_minimal(base_size = 12) +
  labs(
    x = "Sample",
    y = "Number of Resistance Genes",
    fill = "Antibiotic Class",
    title = "Antimicrobial Resistance Genes Across Samples"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    panel.grid.major.x = element_blank()
  )

ggsave(
  filename = "stacked_barplt.svg",
  plot = stacked_bar,
  device = "svg",
  width = 10,
  height = 6
)
#ggsave("stacked_barplt.png", plot = stacked_bar, width = 10, height = 6, dpi = 600)

#sort the df according to antibiotic class
df_long <- arrange(df_long, antibiotic)
df_long$Gene <- factor(df_long$Gene, levels = unique(df_long$Gene))


#head(df_long)
dflong1 <- df_long
colnames(dflong1)[4] <- "Drug_Class"
htmp <- ggplot(dflong1, aes(x = sample, y = Gene)) +
  geom_tile(aes(alpha = Presence,fill=Drug_Class),color="white") +
  scale_alpha(range = c(0, 1), guide = "none") + # fade out absent values
  scale_x_discrete(expand = c(0, 0)) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.5, size = 12, face= "bold"),
    axis.title = element_blank(),
    legend.title = element_text(face = "bold", size = 12),
    legend.text =  element_text(face = "bold", size = 10),
    axis.text.y = element_text(size = 11, face= "bold"),
    #legend.position = "bottom",
    plot.margin = margin(t = 5.5, r = 5.5, b = 60, l = 5.5)
   
  ) 

ggsave(
  filename = "htmap_gene_profiles.svg",
  plot = htmp,
  device = "svg",
  width = 10,
  height = 10
  
)


library(cowplot)
htmp_nolgnd <- htmp + theme(legend.position = "none")
ggsave(
  filename = "htmp_nolgnd.svg",
  plot = htmp_nolgnd,
  device = "svg",
  width = 10,
  height = 10
  
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


#ggsave("htmap_gene_profiles.png", plot = htmp, width = 10, height = 6, dpi = 600)





# sts$strain
# sts$escherichia__mlst_achtman__ST
# sts$escherichia__ezclermont__Clermont_type
# sts$escherichia__ectyper__Serotype

sts_sub <- sts[,c("strain","escherichia__mlst_achtman__ST","escherichia__ezclermont__Clermont_type",
                  "escherichia__ectyper__Serotype")]

colnames(sts_sub) <- c("strain","ST","Phylogroup","Serotype")
head(sts_sub)
sts_sub$Phylogroup %>% table()

ggplot(sts_sub, aes(x = strain, y = ST, color = Phylogroup)) +
  geom_point(size = 3) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Sample",
    y = "Sequence Type (ST)",
    color = "Phylogroup",
    title = "Sequence Types per Sample"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))



sts_sub <- sts_sub %>%
  mutate(
    ST = factor(ST, levels = unique(ST)),
    Phylogroup = factor(Phylogroup, levels = unique(Phylogroup)),
    strain = factor(strain, levels = sts_sub$strain)  # preserve sample order
  )

# Plot


ggplot(sts_sub, aes(x = strain, y = ST, fill = Phylogroup)) +
  geom_tile(color = "white", width = 0.8, height = 0.8) +  # tiles for clarity
  scale_fill_brewer(palette = "Set2") +                     # publication-friendly palette
  theme_minimal(base_size = 12) +
  labs(
    x = "Sample",
    y = "Sequence Type (ST)",
    fill = "Phylogroup",
    title = "ST and Phylogroup Assignment per Sample"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 10, face = "bold"),
    panel.grid = element_blank()
  )









df_long$sample <- as.character(df_long$sample)
sts_sub$strain <- as.character(sts_sub$strain)

df_long2 <- df_long %>%
  left_join(sts_sub, by = c("sample"="strain"))

df_long2 <- df_long2 %>%
  mutate(Phylogroup = factor(Phylogroup, levels = unique(Phylogroup)))

df_long2$sample <- factor(df_long2$sample, levels = ind_inc_sample )

phylo_strip <- df_long2 %>%
  distinct(sample, Phylogroup) %>%
  group_by(sample) %>%
  summarize(ypos = max(df_long2$Presence) + 0.5) %>%  # offset above bars
  left_join(distinct(df_long2, sample, Phylogroup), by = "sample")





p <- ggplot(df_long2, aes(x = sample, y = Presence, fill = antibiotic)) +
  geom_bar(stat = "identity") +
  geom_point(data = phylo_strip, aes(x = sample, y = ypos, shape = Phylogroup),
             size = 3, color = "black", inherit.aes = FALSE) +
  scale_fill_manual(values = pub_palette_18) +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.1, size = 14,face = "bold"),
        axis.text.y = element_text(size = 15,face = "bold"),
        legend.title = element_text(face = "bold", size = 12),
        #legend.position = c(0.5, 0.5),
        #legend.justification = "center",
        #legend.key.size = unit(0.1, "cm"),
        legend.text  = element_text(size = 10),
        title = element_text(face = "bold", size = 20),
        
          
        panel.grid = element_blank()) +
  
  labs(fill = "Antibiotic Class",
       shape = "Phylogroup",
       y = "Number of AMR Genes",
       x = "Sample",
       title = "AMR Gene Profiles with Phylogroup Shapes") 
p <- p +
  scale_x_discrete(expand = c(0, 0)) +
  theme(
    panel.spacing = unit(0, "pt"),
    plot.margin = unit(c(0.5, 0.5, 0, 0.5), "cm")
  )
p_withlegend <- p

ggsave(
  filename = "AMR_gene_profiles_Pylogroups.svg",
  plot =p,
  device = "svg",
  width = 16,
  height = 10
  
)
##plot minus legend
p <- p + theme(legend.position = "none")
ggsave(
  filename = "AMR_gene_profiles_Pylogroups_Nolegend.svg",
  plot =p,
  device = "svg",
  width = 16,
  height = 10
  
)


## getting the legends

library(cowplot)
legend_only <- get_legend(p_withlegend)
ggsave("legend_only_Barwithpylogroups.svg", legend_only, width = 10, height = 10)



#ggsave("AMR_gene_profiles.png", plot = p, width = 10, height = 6, dpi = 600)


table(sts_sub$Serotype)

####Stacked bar plot for antibiotic classes


drugclas_freq <- as.data.frame(drugclas_freq) 
drugclas_freq$Freq %>% sort(decreasing = TRUE)
drugclas_freq <- drugclas_freq[order(drugclas_freq$Freq,decreasing = TRUE), ]
class_freq <- drugclas_freq$.

write.table(file="Antibiotic_Family.txt",drugclas_freq)

df_long2$antibiotic <- factor(df_long2$antibiotic, levels = drugclas_freq$.)



ggplot(df_long2, aes(x = antibiotic,y=Presence,fill=antibiotic)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8,face = "bold"),
        axis.text.y = element_text(size = 10,face = "bold"),
        legend.title = element_text(face = "bold", size = 12),
        panel.grid = element_blank()) +
  labs(
       y = "Number of Genes",
       x = "Antibiotic Class",
       title = "Frequency of the different Antibiotic Classes with Resistance")


df_long2 <- left_join(df_long2,drugclas_freq,join_by(antibiotic==.), relationship ="many-to-one")

ggplot(df_long2, aes(x = antibiotic, y = Presence, fill = antibiotic)) +
  geom_col(show.legend = FALSE) +
  geom_text(
    aes(label = Freq),
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
    title = "Frequency of Antibiotic Classes with Resistance"
  )



df_class <- df_long2 %>%
  group_by(antibiotic) %>%
  summarise(
    Presence = sum(Presence),
    Freq = sum(Freq)
  )

df_long_uniqueG <- df_long2[!duplicated(df_long2$Gene),]
uniqueG <- df_long_uniqueG$antibiotic %>% table() %>% as.data.frame()
colnames(uniqueG)[2] <- "Unique_Gene_Counts"

df_long2 <- left_join(df_long2,uniqueG,join_by(antibiotic==.), relationship ="many-to-one")

df_class <- left_join(df_class,uniqueG,join_by(antibiotic==.), relationship ="many-to-one")

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









noG_sample <- genes_amr$Name %>% table() %>% as.data.frame()
noG_sample <- noG_sample[order(noG_sample$Freq, decreasing = TRUE),]
summary(noG_sample$Freq)



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
colnames(anti_B_p_b )[which(colnames(anti_B_p_b)=="LINCOSAMIDE/MACROLIDE/STREPTOGRAMIN")] <- "LIN/MAC/STRG"
colnames(anti_B_p_b )[which(colnames(anti_B_p_b)=="AMINOGLYCOSIDE/QUINOLONE")] <- "AMG/QUI" 
colnames(anti_B_p_b )[which(colnames(anti_B_p_b)=="PHENICOL/QUINOLONE")] <- "PHE/QUI"
colnames(anti_B_p_b )[which(colnames(anti_B_p_b)=="MACROLIDE/STREPTOGRAMIN")] <- "MAC/STRG"
m <- apply(anti_B_p_b[,-1],2,as.numeric)

apply(m,2,sum) %>% as.data.frame() ##frequencies of the diffent antibiotic classes


ord_antiBclass <- apply(anti_B_p_b[-1],2,sum)  %>% as.data.frame() 

ord_antiBclass <- rownames(ord_antiBclass)[order(ord_antiBclass$.,decreasing = TRUE)]
anti_B_p_b <-  cbind(anti_B_p_b$sample, anti_B_p_b[,ord_antiBclass]) ## arranged according to resistance
colnames(anti_B_p_b)[1] <- "sample"
write.csv(anti_B_p_b,file="anti_B_presence.csv", row.names = FALSE)


write.csv(anti_B_p_b,file="anti_B_presence.csv", row.names = FALSE)
gn <- apply(anti_B_p_b[,-1],2,as.numeric)
anti_B_p_b[,-1] <-  gn



df_long_present <- df_long[df_long$Presence==1,]


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



df_long_present <- df_long_small[df_long_small$Presence==1,]
subclas <- df_long_present$subclass %>% table() %>% as.data.frame()
ord_subclas <- subclas$.[order(subclas$Freq,decreasing = TRUE)]
df_long$subclass <- factor(df_long$subclass, levels = ord_subclas) ##order acording to number subclass

sample_ord <- df_long_small$sample %>% table() %>% as.data.frame() ##order accroding to number sample
ord_sample <- sample_ord$.[order(sample_ord$Freq,decreasing = TRUE)]
df_long$sample <- factor(df_long$sample, levels = ord_sample)



plt_subclass <- ggplot(df_long, aes(x = sample, y = subclass)) +
  geom_tile(aes(alpha = Presence),fill = "red", color="white") +
  scale_alpha(range = c(0, 1), guide = "none") + # fade out absent values
  scale_x_discrete(expand = c(0, 0)) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.5, size = 6, face= "bold"),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 8, face= "bold")
  )+
  labs(title = "Resistance Presence or absence per individual Antibiotic")

ggsave("burden_per_antibiotic_eoli.png", plot = plt_subclass, width = 12, height = 6, dpi = 600)




#write.csv(anti_B_p_b,file="anti_B_presence.csv", row.names = FALSE)














#making  color strip for the phylogroups

phy_num <- unique(sts_sub$Phylogroup)

cols20 <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173",
  "#5254a3", "#9c9ede", "#6b6ecf", "#e7969c", "#31a354"
)

phylogroup_colors <- c(
  "A"  = "#1b9e77",  # green
  "B1" = "#d95f02",  # orange
  "B2" = "#7570b3",  # purple
  "C"  = "#e7298a",  # pink
  "D"  = "#66a61e",  # olive green
  "E"  = "#e6ab02",  # mustard
  "F"  = "#a6761d",  # brown
  "G"  = "#666666",  # gray
  "U" = "#1f77b4"
)

phylogroups <- c("A", "B1", "B2", "C", "D", "E", "F", "G","U")

as.list(phylogroup_colors)
sts_sub2 <- sts_sub

sts_sub <- sts_sub[sts_sub$Phylogroup!="EC_control_fail",]


sts_sub$color_phylo <- NA

for (i in 1:nrow(sts_sub)){
  ind <- (phylogroups==sts_sub$Phylogroup[i]) %>% which()
  #print(ind)
  grp <- phylogroups[ind]
  sts_sub$color_phylo[i] <- phylogroup_colors[[grp]]
}

write.csv(sts_sub[,c(1,5,3)], file = "color_strip_phylo.csv",row.names = FALSE)


### clonal complexes

library(tidyverse)
cplx <- sts$escherichia__mlst_achtman__clonal_complex %>% gsub("ST|Cplx","", x=.) #(pattern, replacement, data)
cplx <- paste("CC",cplx,sep = "")

sts_sub2$CC <- cplx
sts_sub2$CC[sts_sub2$CC=="CC"] <- NA
sts_sub2$CC[sts_sub2$CC=="CC-"] <- NA

cc_vector <- sts_sub2$CC[!is.na(sts_sub2$CC)]

n <- length(unique(cc_vector))
cc_colors <- setNames(scales::hue_pal()(n), unique(cc_vector)) ##set automatic colors
cc_colors <- as.list(cc_colors)
names(cc_colors) <- trimws(names(cc_colors)) ##removes spaces
cc_colors[["CC69"]]

cc_vector <- unique(cc_vector)


sts_sub2$CC_color <- NA
sts_sub2$CC_itol <- NA
for (i in 1:nrow(sts_sub2)){
  if(!is.na(sts_sub2$CC)[i]){
    ind <- (cc_vector==sts_sub2$CC[i]) %>% which()
    cplx <- cc_vector[ind]
    cplx <- trimws(as.character(sts_sub2$CC[i]))
    print(cplx)
    sts_sub2$CC_color[i] <- cc_colors[[cplx]]
    itol <- paste(sts_sub2$strain[i],"label","leaf",cc_colors[[cplx]],"1","normal", sep = ",")
    sts_sub2$CC_itol[i] <- itol
  } else {
    itol <- paste(sts_sub2$strain[i],"label","leaf","#000000","1","normal", sep = ",")
    sts_sub2$CC_itol[i] <- itol}
}
cc_itol <- sts_sub2$CC_itol

#write.csv(cc_itol[!is.na(cc_itol)],file="clonal_complex.csv")
write.csv(cc_itol,file="clonal_complex.csv",row.names = FALSE)


ggplot(sts_sub2, aes(x = strain, y = CC, fill = Phylogroup)) +
  geom_tile(color = "white", width = 0.8, height = 0.8) +  # tiles for clarity
  scale_fill_brewer(palette = "Set2") +                     # publication-friendly palette
  theme_minimal(base_size = 12) +
  labs(
    x = "Sample",
    y = "Sequence Type (ST)",
    fill = "Phylogroup",
    title = "ST and Phylogroup Assignment per Sample"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 10, face = "bold"),
    panel.grid = element_blank()
  )




###PLASMIDS

plasmid1 <- read.csv(file="~/marina/prof_dam/ecoli/plasmids_matrix.txt",
                     sep="\t", header=TRUE)
plasmid1$X.FILE <- gsub(".fasta","",plasmid1$X.FILE )
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

plasmid_only <- sapply(plasmid1[,-c(1,2)], pre_abs)
plasmid1 <- cbind(V1,plasmid_only)



colnames(plasmid1)[-c(1,2)]



#plas_fam <- read.csv2(file="~/Downloads/plasmid_family_classification (1).csv", sep = ",", header = TRUE)

family <- c(
  "Col","Col","Col","Col","Col","Col","Col","Col","Col","Col","Col",
  "IncB/O/K/Z","IncB/O/K/Z","IncB/O/K/Z",
  "IncF","IncF","IncF","IncF","IncF","IncF","IncF","IncF","IncF","IncF","IncF","IncF",
  "IncHI1B","IncI1",
  "IncX","IncX","IncY","p0111"
)





colnames(plasmid1)[-c(1,2)]

# plasmid_df <- data.frame(
#   plasmid = plasmid1,
#   family = family
# )

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


ggsave("Plasmid.png", plot = plasmid_plt, width = 10, height = 10, dpi = 600)

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

sym <- c(1,2,3,4,5,6) ##adjust depending on point antibiotics unique
sym <- list(1,2,3,4,5,6)
names(sym) <- points_ant

col_pts <- list("#ff0000","#00ff00","#0000ff","#cfd623","#00CED1","#FF7F00")


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




amr_gene_freq <- apply(gene_p_b[,-1],1,sum)

paste(gene_p_b$sample,amr_gene_freq )






itol_color_grad <- paste(gene_p_b$sample,amr_gene_freq) %>% as.data.frame()


write.csv(file="colorstrip_number_resistance.csv",itol_color_grad, row.names = FALSE)

itol_color_grad2 <- cbind(gene_p_b$sample,amr_gene_freq) %>% as.data.frame()
itol_color_grad2$amr_range <- NA
itol_color_grad2$amr_gene_freq <- as.numeric(itol_color_grad2$amr_gene_freq)

for (i in 1:nrow(itol_color_grad2)){
  if(itol_color_grad2$amr_gene_freq[i] <= 9){
    itol_color_grad2$amr_range[i] <- 9
  } else if((itol_color_grad2$amr_gene_freq[i]) >9 & (itol_color_grad2$amr_gene_freq[i] < 18)){
    itol_color_grad2$amr_range[i] <- 18
  } else{itol_color_grad2$amr_range[i] <- 27}
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

colors <- c("#FFFF33","#009E73","#654321", "#E69F00","#D55E00", "#CC79A7", "#56B4E9" )


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
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    panel.grid.major.x = element_blank()
  )
ggsave("stacked_bar_points.png", plot = bar_points, width = 15, height = 6, dpi = 600)



### frequencies of point muts


frq_pts <- all_point$Element.symbol %>% table() %>% as.data.frame()

pts_clas <- all_point$Class[!duplicated(all_point$Element.symbol)]
write.table(cbind (frq_pts,pts_clas),file="freq_points.txt")



##number of point mutations
pts <- all_point$Element.symbol %>% table() %>% as.data.frame()
colnames(pts)[1] <-"Element.symbol" 

k <- all_point[!duplicated(all_point$Element.symbol),c(1,7,13)]

k <- merge(pts,k,"Element.symbol")

write.table(k[,c(1,2,4)],file="pts_freq.txt")
