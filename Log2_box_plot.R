setwd("C:/Users/ASUS/OneDrive/Desktop/Chili")
library(GEOquery)
library(stringr)
library(edgeR)
library(tidyverse)

#### load the count_data ####
count_data<- as.data.frame(read.table("GSE132824_Abiotic_RNA-Seq.Readcount.txt"))

#### load the col_data ####
gse<- getGEO(GEO= "GSE132824", GSEMatrix = TRUE)
metadata<- pData(phenoData(gse[[1]]))

col_data<- metadata%>%
  select(1,44,45)%>%
  as.data.frame()%>%
  dplyr::rename(
    stress = `stress:ch1`,
    time = `time:ch1`
  )%>%
  dplyr:: mutate(title=gsub("Leaves,","",title))%>%
  mutate(title=gsub("Abio-","",title))%>%
  mutate(title=gsub("-",".",title))%>%
  mutate(title= str_to_lower(title))%>%
  mutate(title= gsub("nacl","NaCl", title))%>%
  mutate(stress=gsub("Abio-mock","Abio_mock",stress))%>%
  mutate(title=str_trim(title, side = "left"))

rownames(col_data)<- NULL
names(col_data)[1]<- "Samples"

# save the col_data
write_csv(col_data, file = "metadata.csv")

# arrange the count_data according to col_data
count_data<- count_data%>%
  select(all_of(col_data$Samples))

# observe the data 
log2_df<- count_data%>%
  gather(key = "Samples", value = "Counts")%>%
  mutate(Samples= factor(Samples, levels = col_data$Samples))%>%
  mutate(Sample_type =sub("\\..*$","",Samples))%>%
  mutate(Sample_type= factor(Sample_type, levels = c("mock","cold","heat","man","NaCl")))%>%
  mutate(log_Counts= log2(Counts+1))
  
# plot the boxplot
ggplot(log2_df, aes(Samples, log_Counts, fill = Sample_type)) +
  geom_boxplot() +
  facet_grid(~ Sample_type, scales = "free_x", space = "free_y", switch = "x") +
  scale_fill_manual(values = c("darkgray","steelblue","red","gold","green3")) +
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    strip.placement = "outside",
    strip.text.x = element_text(size = 12, face = "bold"),
    legend.position = "none",         
    axis.title.x = element_blank(),
    panel.spacing = unit(0, "cm"),
    panel.background = element_blank(),
    plot.background = element_rect(colour = "black", fill = NA, linewidth = 1), # border around whole plot
    strip.background = element_rect(fill = "white", colour = "black",linewidth = 0.8)) +
  ylab("Log2 read numbers")

