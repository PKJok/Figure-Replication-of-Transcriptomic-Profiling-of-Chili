
setwd("C:/Users/ASUS/OneDrive/Desktop/Chili")

library(stringr)
library(edgeR)
library(tidyverse)
library(limma)
library(gridExtra)

# load the count_data
count_data<- as.data.frame(read.table("GSE132824_Abiotic_RNA-Seq.Readcount.txt"))

col_data<- read.csv("metadata.csv")
col_data<- col_data%>%
  column_to_rownames(var = "Samples")
dim(count_data)


count_data<- count_data%>%
  select(all_of(rownames(col_data)))

############### Preprocessing ###################
# creating different data frames of stresses
cold_df<- count_data%>%
  select(1:33)

heat_df<- count_data%>%
  select(1:18,34:48)

NaCl_df<- count_data%>%
  select(1:18,49:63)

man_df<-count_data%>%
  select(1:18,64:78)

# Select the coldata for each stress
cold_colData<- col_data[c(1:33),]
heat_colData<- col_data[c(1:18,34:48),]
NaCl_colData<- col_data[c(1:18,49:63),]
man_colData<- col_data[c(1:18,64:78),]


#### Shapes for Timepoint ####
my_shapes<- c(1,2,0,16,17,15)

 

#----------------------------------------------------------------------#

####------------------------ Heat ---------------------------------#####

#----------------------------------------------------------------------#

y_heat<- DGEList(heat_df)
dim(y_heat)

# adding group to the sample
heat_group<- paste(heat_colData$stress,heat_colData$time, sep = ".")
group<- factor(heat_group)

y_heat$samples$group<- group

# filter the lowly expressed genes
heat_keep<- filterByExpr(y_heat, group = group)
y_heat<- y_heat[heat_keep, ,keep.lib.size= FALSE]
dim(y_heat)

# TMM Normalization 
heat_TMM<- calcNormFactors(y_heat,method = "TMM")
heat_logcount<- cpm(heat_TMM, log = TRUE)

##### "HEAT PCA" #####
heat_pca<- prcomp(t(heat_logcount))
plot(heat_pca$x[,1], heat_pca$x[,2])

heat_pca_per<- round((heat_pca$sdev^2/sum(heat_pca$sdev^2))*100,2)
heat_pca_per[2]


p1<-heat_pca$x %>% 
  as.data.frame() %>%
  select(PC1, PC2) %>%
  rownames_to_column(var = "Samples") %>%
  mutate(
    Treatment = case_when(
      str_detect(Samples, "mock") ~ "mock",
      str_detect(Samples, "cold") ~ "cold",
      str_detect(Samples, "heat") ~ "heat",
      str_detect(Samples, "man") ~ "man",
      str_detect(Samples, "NaCl") ~ "NaCl",
      TRUE ~ "unknown"  # Handle unexpected cases
    ),
    Timepoint = case_when(
      str_detect(Samples, "0h") ~ "0h",
      str_detect(Samples, "3h") ~ "3h",
      str_detect(Samples, "6h") ~ "6h", 
      str_detect(Samples, "12h") ~ "12h",
      str_detect(Samples, "24h") ~ "24h",
      str_detect(Samples, "72h") ~ "72h",
      TRUE ~ "unknown"
    )
  )%>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("0h","3h","6h","12h","24h","72h")),
    Treatment = factor(Treatment, levels = c("mock", "heat", "cold", "man", "NaCl"))
  )%>%
  ggplot(aes(PC1, PC2, colour = Treatment, shape = Timepoint)) +
  geom_point(size = 3, stroke=1.2) +
  #stat_ellipse(aes(color = Treatment), level = 0.95, linewidth = 1) +
  scale_shape_manual(values = my_shapes) +
  scale_color_manual(values = c('mock' = 'gray48',
                                'man' = '#feb81a',
                                'cold' = 'royalblue', 
                                'heat' = 'red',
                                'NaCl' = '#65a947')) +
  theme_bw() +
  labs(x = paste0("PC1 (",heat_pca_per[1], "%)"),
       y = paste0("PC2 (",heat_pca_per[2], "%)")) +
  coord_cartesian(ylim = c(-100, 100), xlim = c(-100, 100)) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(face = "bold",hjust = 0.5))+
  ggtitle("Heat")

#----------------------------------------------------------------------#

####------------------------ Cold ---------------------------------#####

#----------------------------------------------------------------------#

y_cold<- DGEList(cold_df)
dim(y_cold)
y_cold$samples

# adding group info
cold_group<- paste(cold_colData$stress, cold_colData$time, sep = ".")
group<- factor(cold_group)
group

y_cold$samples$group<- group
y_cold$samples

#filter the lowly-expressed genes
cold_keep<- filterByExpr(y_cold, group = group)
y_cold<- y_cold[cold_keep, , keep.lib.sizes=FALSE]
dim(y_cold)

# TMM Normalization
cold_TMM<- calcNormFactors(y_cold, method ="TMM")
cold_logcount<- cpm(cold_TMM, log = TRUE)


####### "COLD PCA" #######
cold_pca<- prcomp(t(cold_logcount))
plot(cold_pca$x[,1], cold_pca$x[,2])

cold_pca_per<- round((cold_pca$sdev^2/sum(cold_pca$sdev^2))*100,2)
  
p2<- cold_pca$x %>% 
  as.data.frame() %>%
  select(PC1, PC2) %>%
  rownames_to_column(var = "Samples") %>%
  mutate(
    Treatment = case_when(
      str_detect(Samples, "mock") ~ "mock",
      str_detect(Samples, "cold") ~ "cold",
      str_detect(Samples, "heat") ~ "heat",
      str_detect(Samples, "man") ~ "man",
      str_detect(Samples, "NaCl") ~ "NaCl",
      TRUE ~ "unknown"  # Handle unexpected cases
    ),
    Timepoint = case_when(
      str_detect(Samples, "0h") ~ "0h",
      str_detect(Samples, "3h") ~ "3h",
      str_detect(Samples, "6h") ~ "6h", 
      str_detect(Samples, "12h") ~ "12h",
      str_detect(Samples, "24h") ~ "24h",
      str_detect(Samples, "72h") ~ "72h",
      TRUE ~ "unknown"
    )
  ) %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("0h","3h","6h","12h","24h","72h")),
    Treatment = factor(Treatment, levels = c("mock", "heat", "cold", "man", "NaCl"))
  )%>%
  ggplot(aes(PC1, PC2, colour = Treatment, shape = Timepoint)) +
  geom_point(size = 3, stroke=1.2) +
  scale_shape_manual(values = my_shapes) +
  scale_color_manual(values = c('mock' = 'gray48',
                                'man' = '#feb81a',
                                'cold' = 'royalblue', 
                                'heat' = 'red',
                                'NaCl' = '#65a947')) +
  theme_bw() +
  labs(x = paste0("PC1 (",cold_pca_per[1], "%)"),
       y = paste0("PC2 (",cold_pca_per[2], "%)")) +
  coord_cartesian(ylim = c(-100, 100), xlim = c(-100, 150)) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(face = "bold",hjust = 0.5))+
  ggtitle("Cold")

#----------------------------------------------------------------------#

####------------------------ Mannitol -----------------------------#####

#----------------------------------------------------------------------#

y_man<- DGEList(man_df)
dim(y_man)
y_man$samples

# adding group info
man_group<- paste(man_colData$stress, man_colData$time, sep = ".")
group<- factor(man_group)
group

y_man$samples$group<- group
y_man$samples

#filter the lowly-expressed genes
man_keep<- filterByExpr(y_man, group = group)
y_man<- y_man[man_keep, , keep.lib.sizes=FALSE]
dim(y_man)

# TMM Normalization
man_TMM<- calcNormFactors(y_man, method ="TMM")
man_logcount<- cpm(man_TMM, log = TRUE)


##### "MANNITOL PCA" ####

man_pca<- prcomp(t(man_logcount))

man_pca_per<- round((man_pca$sdev^2/sum(man_pca$sdev^2))*100,2)

p3<- man_pca$x %>% 
  as.data.frame() %>%
  select(PC1, PC2) %>%
  rownames_to_column(var = "Samples") %>%
  mutate(
    Treatment = case_when(
      str_detect(Samples, "mock") ~ "mock",
      str_detect(Samples, "cold") ~ "cold",
      str_detect(Samples, "heat") ~ "heat",
      str_detect(Samples, "man") ~ "man",
      str_detect(Samples, "NaCl") ~ "NaCl",
      TRUE ~ "unknown"  # Handle unexpected cases
    ),
    Timepoint = case_when(
      str_detect(Samples, "0h") ~ "0h",
      str_detect(Samples, "3h") ~ "3h",
      str_detect(Samples, "6h") ~ "6h", 
      str_detect(Samples, "12h") ~ "12h",
      str_detect(Samples, "24h") ~ "24h",
      str_detect(Samples, "72h") ~ "72h",
      TRUE ~ "unknown"
    )
  ) %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("0h","3h","6h","12h","24h","72h")),
    Treatment = factor(Treatment, levels = c("mock", "heat", "cold", "man", "NaCl"))
  )%>%
  ggplot(aes(PC1, PC2, colour = Treatment, shape = Timepoint)) +
  geom_point(size = 3, stroke=1.2) +
  scale_shape_manual(values = my_shapes) +
  scale_color_manual(values = c('mock' = 'gray48',
                                'man' = '#feb81a',
                                'cold' = 'royalblue', 
                                'heat' = 'red',
                                'NaCl' = '#65a947')) +
  theme_bw() +
  labs(x = paste0("PC1 (", man_pca_per[1], "%)"),
       y = paste0("PC2 (", man_pca_per[2], "%)")) +
  coord_cartesian(ylim = c(-100, 100), xlim = c(-100, 100)) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(face = "bold",hjust = 0.5))+
  ggtitle("Man")


#----------------------------------------------------------------------#

####-------------------------- NaCl -------------------------------#####

#----------------------------------------------------------------------#

y_NaCl<- DGEList(NaCl_df)
dim(y_NaCl)

# adding group to the y_NaCl
NaCl_group<- paste(NaCl_colData$stress, NaCl_colData$time, sep = ".")
group<- factor(NaCl_group)

y_NaCl$samples$group<- group

# filtering lowly expressed genes
NaCl_keep<- filterByExpr(y_NaCl, group = group)
y_NaCl<- y_NaCl[NaCl_keep, keep.lib.sizes= FALSE]
dim(y_NaCl)

# TMM Normalization
NaCl_TMM<- calcNormFactors(y_NaCl, method = "TMM")
NaCl_logcount<- cpm(NaCl_TMM, log = TRUE)

##### "NaCl PCA" ####
NaCl_pca<- prcomp(t(NaCl_logcount))

NaCl_pca_per<- round((NaCl_pca$sdev^2/sum(NaCl_pca$sdev^2))*100,2)

p4<- NaCl_pca$x %>% 
  as.data.frame() %>%
  select(PC1, PC2) %>%
  rownames_to_column(var = "Samples") %>%
  mutate(
    Treatment = case_when(
      str_detect(Samples, "mock") ~ "mock",
      str_detect(Samples, "cold") ~ "cold",
      str_detect(Samples, "heat") ~ "heat",
      str_detect(Samples, "man") ~ "man",
      str_detect(Samples, "NaCl") ~ "NaCl",
      TRUE ~ "unknown"  # Handle unexpected cases
    ),
    Timepoint = case_when(
      str_detect(Samples, "0h") ~ "0h",
      str_detect(Samples, "3h") ~ "3h",
      str_detect(Samples, "6h") ~ "6h", 
      str_detect(Samples, "12h") ~ "12h",
      str_detect(Samples, "24h") ~ "24h",
      str_detect(Samples, "72h") ~ "72h",
      TRUE ~ "unknown"
    )
  ) %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("0h","3h","6h","12h","24h","72h")),
    Treatment = factor(Treatment, levels = c("mock", "heat", "cold", "man", "NaCl"))
  )%>%
  ggplot(aes(PC1, PC2, colour = Treatment, shape = Timepoint)) +
  geom_point(size = 3, stroke=1.2) +
  scale_shape_manual(values = my_shapes) +
  scale_color_manual(values = c('mock' = 'gray48',
                                'man' = '#feb81a',
                                'cold' = 'royalblue', 
                                'heat' = 'red',
                                'NaCl' = '#65a947')) +
  theme_bw() +
  labs(x = paste0("PC1 (",NaCl_pca_per[1], "%)"),
       y = paste0("PC2 (", NaCl_pca_per[2], "%)")) +
  coord_cartesian(ylim = c(-100, 100), xlim = c(-100, 100)) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(face = "bold",hjust = 0.5))+
  ggtitle("NaCl")


################ plot the combined plot ##############
grid.arrange(p1,p2,p3,p4, nrow=2, ncol=2)












