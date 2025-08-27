# improved MD Plots
setwd("C:/Users/ASUS/OneDrive/Desktop/Chili")

# load the libraries
library(readxl)
library(tidyverse)
library(DESeq2)
library(edgeR)
library(limma)

# load the files
list.files("Normalized TMM")
list.files("C:/Users/ASUS/OneDrive/Desktop/Chili/DEG result for each stress")

########## load DEG of each stress #######

# DEG of heat stress
DEG_heat<- read_excel("DEG result for each stress/DEG result with heat stress.xlsx",
                      skip = 1)%>%
  select(-c(1,3:5))

colnames(DEG_heat)<- DEG_heat[1,]
DEG_heat<- DEG_heat[-1,]

# DEG of Cold Stress
DEG_cold<- read_excel("DEG result for each stress/DEG result with cold stress.xlsx",
                      skip = 1)%>%
  select(-c(1,3:5))

colnames(DEG_cold)<- DEG_cold[1,]
DEG_cold<- DEG_cold[-1,]

# DEG of Mannitol stress
DEG_man<- read_excel("DEG result for each stress/DEG result with man stress.xlsx",
                      skip = 1)%>%
  select(-c(1,3:5))

colnames(DEG_man)<- DEG_man[1,]
DEG_man<- DEG_man[-1,]

# DEG of NaCl Stress
DEG_NaCl<- read_excel("DEG result for each stress/DEG result with NaCl stress.xlsx",
                     skip = 1)%>%
  select(-c(1,3:5))

colnames(DEG_NaCl)<- DEG_NaCl[1,]
DEG_NaCl<- DEG_NaCl[-1,]


###### load the TMM normalized data of heat stress #####
# TMM of Heat
TMM_countData_heat<- read_excel("Normalized TMM/Expression profiling of TMM normalization on heat stress.xlsx")
names(TMM_countData_heat)[1]<- "Gene_ID"

TMM_countData_heat<- TMM_countData_heat%>%
  mutate(Gene_ID=gsub("TC.","",Gene_ID))

# TMM of cold
list.files("Normalized TMM") # extra space in the file name
TMM_countData_cold<- read_excel("Normalized TMM/Expression profiling of TMM normalization on  cold stress.xlsx")
names(TMM_countData_cold)[1]<- "Gene_ID"

TMM_countData_cold<- TMM_countData_cold%>%
  mutate(Gene_ID=gsub("TC.","",Gene_ID))

# TMM of Mannitol
TMM_countData_man<- read_excel("Normalized TMM/Expression profiling of TMM normalization on osmotic stress.xlsx")
names(TMM_countData_man)[1]<- "Gene_ID"

TMM_countData_man<- TMM_countData_man%>%
  mutate(Gene_ID=gsub("TC.","",Gene_ID))

# TMM of NaCl
TMM_countData_NaCl<- read_excel("Normalized TMM/Expression profiling of TMM normalization on NaCl stress.xlsx")
names(TMM_countData_NaCl)[1]<- "Gene_ID"

TMM_countData_NaCl<- TMM_countData_NaCl%>%
  mutate(Gene_ID=gsub("TC.","",Gene_ID))



####### Pre-procesing #######
# DEG processing
# subset for time points of heat stress 
DEG_heat_3h<- DEG_heat[,c(1:3)]
colnames(DEG_heat_3h)[1]<- "Gene_ID"
DEG_heat_3h$Time<- "3h"

DEG_heat_6h <- DEG_heat[,c(1,5,6)]
colnames(DEG_heat_6h)[1]<- "Gene_ID"
DEG_heat_6h$Time<- "6h"

DEG_heat_12h <- DEG_heat[,c(1,8,9)]
colnames(DEG_heat_12h)[1]<- "Gene_ID"
DEG_heat_12h$Time <- "12h"

DEG_heat_24h <- DEG_heat[,c(1,11,12)]
colnames(DEG_heat_24h)[1]<- "Gene_ID"
DEG_heat_24h$Time <- "24h"

DEG_heat_72h <- DEG_heat[,c(1,14,15)]
colnames(DEG_heat_72h)[1]<- "Gene_ID"
DEG_heat_72h$Time <- "72h"


# subset for time points of cold stress 
DEG_cold_3h<- DEG_cold[,c(1:3)]
colnames(DEG_cold_3h)[1]<- "Gene_ID"
DEG_cold_3h$Time<-"3h"

DEG_cold_6h <- DEG_cold[,c(1,5,6)]
colnames(DEG_cold_6h)[1]<- "Gene_ID"
DEG_cold_6h$Time <- "6h"

DEG_cold_12h <- DEG_cold[,c(1,8,9)]
colnames(DEG_cold_12h)[1]<- "Gene_ID"
DEG_cold_12h$Time <- "12h"

DEG_cold_24h <- DEG_cold[,c(1,11,12)]
colnames(DEG_cold_24h)[1]<- "Gene_ID"
DEG_cold_24h$Time <- "24h"

DEG_cold_72h <- DEG_cold[,c(1,14,15)]
colnames(DEG_cold_72h)[1]<- "Gene_ID"
DEG_cold_72h$Time <- "72h"

# subset for time points of man stress 
DEG_man_3h<- DEG_man[,c(1:3)]
colnames(DEG_man_3h)[1]<- "Gene_ID"
DEG_man_3h$Time <- "3h"

DEG_man_6h <- DEG_man[,c(1,5,6)]
colnames(DEG_man_6h)[1]<- "Gene_ID"
DEG_man_6h$Time <- "6h"

DEG_man_12h <- DEG_man[,c(1,8,9)]
colnames(DEG_man_12h)[1]<- "Gene_ID"
DEG_man_12h$Time <- "12h"

DEG_man_24h <- DEG_man[,c(1,11,12)]
colnames(DEG_man_24h)[1]<- "Gene_ID"
DEG_man_24h$Time <- "24h"

DEG_man_72h <- DEG_man[,c(1,14,15)]
colnames(DEG_man_72h)[1]<- "Gene_ID"
DEG_man_72h$Time <- "72h"

# subset for time points of NaCl stress 
DEG_NaCl_3h<- DEG_NaCl[,c(1:3)]
colnames(DEG_NaCl_3h)[1]<- "Gene_ID"
DEG_NaCl_3h$Time <- "3h"

DEG_NaCl_6h <- DEG_NaCl[,c(1,5,6)]
colnames(DEG_NaCl_6h)[1]<- "Gene_ID"
DEG_NaCl_6h$Time <- "6h"

DEG_NaCl_12h <- DEG_NaCl[,c(1,8,9)]
colnames(DEG_NaCl_12h)[1]<- "Gene_ID"
DEG_NaCl_12h$Time <- "12h"

DEG_NaCl_24h <- DEG_NaCl[,c(1,11,12)]
colnames(DEG_NaCl_24h)[1]<- "Gene_ID"
DEG_NaCl_24h$Time <- "24h"

DEG_NaCl_72h <- DEG_NaCl[,c(1,14,15)]
colnames(DEG_NaCl_72h)[1]<- "Gene_ID"
DEG_NaCl_72h$Time <- "72h"


# Process the TMM of each time points 

# TMM of Heat at 3h 
TMM_heat_3h<- TMM_countData_heat[,c(1,20:22)]%>%
  column_to_rownames(var = "Gene_ID")%>%
  mutate(mean_normCount= apply(.,1,mean))%>%
  select(4)%>%
  rownames_to_column(var = "Gene_ID")

# TMM of Heat at 6h 
TMM_heat_6h<- TMM_countData_heat[,c(1,23:25)]%>%
  column_to_rownames(var = "Gene_ID")%>%
  mutate(mean_normCount= apply(.,1,mean))%>%
  select(4)%>%
  rownames_to_column(var = "Gene_ID")

# TMM of Heat at 12h 
TMM_heat_12h<- TMM_countData_heat[,c(1,26:28)]%>%
  column_to_rownames(var = "Gene_ID")%>%
  mutate(mean_normCount= apply(.,1,mean))%>%
  select(4)%>%
  rownames_to_column(var = "Gene_ID")
  
# TMM of Heat at 24h 
TMM_heat_24h<- TMM_countData_heat[,c(1,29:31)]%>%
  column_to_rownames(var = "Gene_ID")%>%
  mutate(mean_normCount= apply(.,1,mean))%>%
  select(4)%>%
  rownames_to_column(var = "Gene_ID")

# TMM of Heat at 72h 
TMM_heat_72h<- TMM_countData_heat[,c(1,32:34)]%>%
  column_to_rownames(var = "Gene_ID")%>%
  mutate(mean_normCount= apply(.,1,mean))%>%
  select(4)%>%
  rownames_to_column(var = "Gene_ID")

# TMM Cold 
# TMM of cold at 3h 
TMM_cold_3h<- TMM_countData_cold[,c(1,20:22)]%>%
  column_to_rownames(var = "Gene_ID")%>%
  mutate(mean_normCount= apply(.,1,mean))%>%
  select(4)%>%
  rownames_to_column(var = "Gene_ID")

# TMM of cold at 6h 
TMM_cold_6h<- TMM_countData_cold[,c(1,23:25)]%>%
  column_to_rownames(var = "Gene_ID")%>%
  mutate(mean_normCount= apply(.,1,mean))%>%
  select(4)%>%
  rownames_to_column(var = "Gene_ID")

# TMM of cold at 12h 
TMM_cold_12h<- TMM_countData_cold[,c(1,26:28)]%>%
  column_to_rownames(var = "Gene_ID")%>%
  mutate(mean_normCount= apply(.,1,mean))%>%
  select(4)%>%
  rownames_to_column(var = "Gene_ID")

# TMM of cold at 24h 
TMM_cold_24h<- TMM_countData_cold[,c(1,29:31)]%>%
  column_to_rownames(var = "Gene_ID")%>%
  mutate(mean_normCount= apply(.,1,mean))%>%
  select(4)%>%
  rownames_to_column(var = "Gene_ID")

# TMM of cold at 72h 
TMM_cold_72h<- TMM_countData_cold[,c(1,32:34)]%>%
  column_to_rownames(var = "Gene_ID")%>%
  mutate(mean_normCount= apply(.,1,mean))%>%
  select(4)%>%
  rownames_to_column(var = "Gene_ID")

# TMM Mannitol
# TMM of man at 3h 
TMM_man_3h<- TMM_countData_man[,c(1,20:22)]%>%
  column_to_rownames(var = "Gene_ID")%>%
  mutate(mean_normCount= apply(.,1,mean))%>%
  select(4)%>%
  rownames_to_column(var = "Gene_ID")

# TMM of man at 6h 
TMM_man_6h<- TMM_countData_man[,c(1,23:25)]%>%
  column_to_rownames(var = "Gene_ID")%>%
  mutate(mean_normCount= apply(.,1,mean))%>%
  select(4)%>%
  rownames_to_column(var = "Gene_ID")

# TMM of man at 12h 
TMM_man_12h<- TMM_countData_man[,c(1,26:28)]%>%
  column_to_rownames(var = "Gene_ID")%>%
  mutate(mean_normCount= apply(.,1,mean))%>%
  select(4)%>%
  rownames_to_column(var = "Gene_ID")

# TMM of man at 24h 
TMM_man_24h<- TMM_countData_man[,c(1,29:31)]%>%
  column_to_rownames(var = "Gene_ID")%>%
  mutate(mean_normCount= apply(.,1,mean))%>%
  select(4)%>%
  rownames_to_column(var = "Gene_ID")

# TMM of man at 72h 
TMM_man_72h<- TMM_countData_man[,c(1,32:34)]%>%
  column_to_rownames(var = "Gene_ID")%>%
  mutate(mean_normCount= apply(.,1,mean))%>%
  select(4)%>%
  rownames_to_column(var = "Gene_ID")

# TMM of NaCl
# TMM of man at 3h 
TMM_NaCl_3h<- TMM_countData_NaCl[,c(1,20:22)]%>%
  column_to_rownames(var = "Gene_ID")%>%
  mutate(mean_normCount= apply(.,1,mean))%>%
  select(4)%>%
  rownames_to_column(var = "Gene_ID")

# TMM of NaCl at 6h 
TMM_NaCl_6h<- TMM_countData_NaCl[,c(1,23:25)]%>%
  column_to_rownames(var = "Gene_ID")%>%
  mutate(mean_normCount= apply(.,1,mean))%>%
  select(4)%>%
  rownames_to_column(var = "Gene_ID")

# TMM of NaCl at 12h 
TMM_NaCl_12h<- TMM_countData_NaCl[,c(1,26:28)]%>%
  column_to_rownames(var = "Gene_ID")%>%
  mutate(mean_normCount= apply(.,1,mean))%>%
  select(4)%>%
  rownames_to_column(var = "Gene_ID")

# TMM of NaCl at 24h 
TMM_NaCl_24h<- TMM_countData_NaCl[,c(1,29:31)]%>%
  column_to_rownames(var = "Gene_ID")%>%
  mutate(mean_normCount= apply(.,1,mean))%>%
  select(4)%>%
  rownames_to_column(var = "Gene_ID")

# TMM of NaCl at 72h 
TMM_NaCl_72h<- TMM_countData_NaCl[,c(1,32:34)]%>%
  column_to_rownames(var = "Gene_ID")%>%
  mutate(mean_normCount= apply(.,1,mean))%>%
  select(4)%>%
  rownames_to_column(var = "Gene_ID")


######## Merging DEG and mean_normCounts #######

# merge heat 
# heat at 3h (mean_norm_counts and log2fold and pvalue)
merged_heat_3h<- merge(DEG_heat_3h,TMM_heat_3h)%>%
  mutate(logCPM= log2(mean_normCount+1))

# heat at 6h (mean_norm_counts and log2fold and pvalue)
merged_heat_6h<- merge(DEG_heat_6h,TMM_heat_6h)%>%
  mutate(logCPM= log2(mean_normCount+1))

# heat at 12h (mean_norm_counts and log2fold and pvalue)
merged_heat_12h<- merge(DEG_heat_12h,TMM_heat_12h)%>%
  mutate(logCPM= log2(mean_normCount+1))

# heat at 24h (mean_norm_counts and log2fold and pvalue)
merged_heat_24h<- merge(DEG_heat_24h,TMM_heat_24h)%>%
  mutate(logCPM= log2(mean_normCount+1))

# heat at 72h (mean_norm_counts and log2fold and pvalue)
merged_heat_72h<- merge(DEG_heat_72h,TMM_heat_72h)%>%
  mutate(logCPM= log2(mean_normCount+1))

# merge cold
# cold at 3h (mean_norm_counts and log2fold and pvalue)
merged_cold_3h<- merge(DEG_cold_3h,TMM_cold_3h)%>%
  mutate(logCPM= log2(mean_normCount+1))

# cold at 6h (mean_norm_counts and log2fold and pvalue)
merged_cold_6h<- merge(DEG_cold_6h,TMM_cold_6h)%>%
  mutate(logCPM= log2(mean_normCount+1))

# cold at 12h (mean_norm_counts and log2fold and pvalue)
merged_cold_12h<- merge(DEG_cold_12h,TMM_cold_12h)%>%
  mutate(logCPM= log2(mean_normCount+1))

# cold at 24h (mean_norm_counts and log2fold and pvalue)
merged_cold_24h<- merge(DEG_cold_24h,TMM_cold_24h)%>%
  mutate(logCPM= log2(mean_normCount+1))

# cold at 72h (mean_norm_counts and log2fold and pvalue)
merged_cold_72h<- merge(DEG_cold_72h,TMM_cold_72h)%>%
  mutate(logCPM= log2(mean_normCount+1))

# Mannitol merge
# man at 3h (mean_norm_counts and log2fold and pvalue)
merged_man_3h<- merge(DEG_man_3h,TMM_man_3h)%>%
  mutate(logCPM= log2(mean_normCount+1))

# man at 6h (mean_norm_counts and log2fold and pvalue)
merged_man_6h<- merge(DEG_man_6h,TMM_man_6h)%>%
  mutate(logCPM= log2(mean_normCount+1))

# man at 12h (mean_norm_counts and log2fold and pvalue)
merged_man_12h<- merge(DEG_man_12h,TMM_man_12h)%>%
  mutate(logCPM= log2(mean_normCount+1))

# man at 24h (mean_norm_counts and log2fold and pvalue)
merged_man_24h<- merge(DEG_man_24h,TMM_man_24h)%>%
  mutate(logCPM= log2(mean_normCount+1))

# man at 72h (mean_norm_counts and log2fold and pvalue)
merged_man_72h<- merge(DEG_man_72h,TMM_man_72h)%>%
  mutate(logCPM= log2(mean_normCount+1))


# merge NaCl
# NaCl at 3h (mean_norm_counts and log2fold and pvalue)
merged_NaCl_3h<- merge(DEG_NaCl_3h,TMM_NaCl_3h)%>%
  mutate(logCPM= log2(mean_normCount+1))

# NaCl at 6h (mean_norm_counts and log2fold and pvalue)
merged_NaCl_6h<- merge(DEG_NaCl_6h,TMM_NaCl_6h)%>%
  mutate(logCPM= log2(mean_normCount+1))

# NaCl at 12h (mean_norm_counts and log2fold and pvalue)
merged_NaCl_12h<- merge(DEG_NaCl_12h,TMM_NaCl_12h)%>%
  mutate(logCPM= log2(mean_normCount+1))

# NaCl at 24h (mean_norm_counts and log2fold and pvalue)
merged_NaCl_24h<- merge(DEG_NaCl_24h,TMM_NaCl_24h)%>%
  mutate(logCPM= log2(mean_normCount+1))

# NaCl at 72h (mean_norm_counts and log2fold and pvalue)
merged_NaCl_72h<- merge(DEG_NaCl_72h,TMM_NaCl_72h)%>%
  mutate(logCPM= log2(mean_normCount+1))


#### Row binding each stress  ####

rbind_heat<- rbind(merged_heat_3h, merged_heat_6h, merged_heat_12h, merged_heat_24h, merged_heat_72h)%>%
  mutate(Stress= "Heat")

rbind_cold<- rbind(merged_cold_3h, merged_cold_6h, merged_cold_12h, merged_cold_24h, merged_cold_72h)%>%
  mutate(Stress= "Cold")

rbind_man<- rbind(merged_man_3h, merged_man_6h, merged_man_12h, merged_man_24h, merged_man_72h)%>%
  mutate(Stress= "Man")

rbind_NaCl<- rbind(merged_NaCl_3h, merged_NaCl_6h, merged_NaCl_12h, merged_NaCl_24h, merged_NaCl_72h)%>%
  mutate(Stress= "NaCl")

#### Now Final Binding ####

final_df<- rbind(rbind_heat, rbind_cold, rbind_man, rbind_NaCl)
names(final_df)[3]<-"pvalue"


final_df$logFC<- as.numeric(final_df$logFC)
final_df$pvalue<- as.numeric(final_df$pvalue)
final_df$Time<- factor(final_df$Time, 
                       levels = c("3h","6h","12h","24h","72h"))
final_df$Stress<-factor(final_df$Stress, 
                        levels = c("Heat","Cold","Man","NaCl"))
final_df$logCPM<- as.numeric(final_df$logCPM)

# adding significance column for 
final_df<- final_df%>%
  mutate(Significance= case_when(
    pvalue < 0.05 & logFC>0 ~"Up",
    pvalue < 0.05 & logFC<0 ~"Down",
    pvalue > 0.05           ~"non-DE"
  ))



#### plotting the final_data ####

ggplot(final_df, aes(logCPM, logFC, colour = Significance)) +
  geom_point(size=0.9) +
  scale_color_manual(values = c("Up" = "red",
                                "Down" = "blue",
                                "non-DE" = "gray2")) +
  facet_grid(Stress ~ Time, scales = "free") +
  theme_test() +
  theme(
    strip.placement = "outside",
    strip.text = element_text(size = 11),
    strip.background = element_blank(),
    panel.background = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  guides(colour = guide_legend(override.aes = list(size = 2.5)))+
  geom_hline(yintercept = c(1, -1), col = "green")+
  labs(
    x = "Mean of normalized counts",
    y = "Log2 fold change",
  )










