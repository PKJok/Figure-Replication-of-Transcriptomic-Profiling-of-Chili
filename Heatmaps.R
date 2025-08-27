getwd()
setwd("C:/Users/ASUS/OneDrive/Desktop/Chili")

# load the libraries
library(pheatmap)
library(tidyverse)
library(edgeR)
library(gridExtra)

list.files("C:/Users/ASUS/OneDrive/Desktop/Chili")

# load the countData and colData
countData<- as.data.frame(read.table("GSE132824_Abiotic_RNA-Seq.Readcount.txt",
                       sep = "\t"))


colData<- read.csv("C:/Users/ASUS/OneDrive/Desktop/Chili/metadata.csv")%>%
  column_to_rownames(var = "Samples")

class(colData)

countData<- countData%>%
  select(all_of(rownames(colData)))

############### Preprocessing ###################

# removing "TC." from every rows of countData to match the DEG list name
rownames(countData)<- gsub("TC.CA.PGAv.1.6.","", rownames(countData))


# creating different data frames of stresses
cold_df<- countData[,1:33]
heat_df<- countData[,c(1:18,34:48)]
NaCl_df<- countData[,c(1:18,49:63)]
man_df<-countData[,c(1:18,64:78)]

# Select the coldata for each stress
cold_colData<- colData[c(1:33),]
heat_colData<- colData[c(1:18,34:48),]
NaCl_colData<- colData[c(1:18,49:63),]
man_colData<- colData[c(1:18,64:78),]


#----------------------------------------------------------------------#

####-------------------------- Heat -------------------------------#####

#----------------------------------------------------------------------#

# Steps for normalization
# step1: form DGEList
y_heat<- DGEList(heat_df)
y_heat$samples

#step2:  forming group 
heat_colData$group <- with(heat_colData, paste(stress,time,sep="."))
heat_colData$group<- factor(heat_colData$group)

# step3: filter the lowly expressed genes
heat_keep<- filterByExpr(y_heat, group = heat_colData$group)
y_heat<- y_heat[heat_keep, ,keep.lib.size= FALSE]
dim(y_heat)

# step4: normalization by TMM method
heat_TMM<- calcNormFactors(y_heat, 
                           group= heat_colData$group,
                           method = "TMM")
heat_logcount<- cpm(heat_TMM, log = TRUE)


# load the DEG data
# DEG of Heat stress
library(readxl)
DEG_heat<- read_excel("DEG result for each stress/DEG result with heat stress.xlsx", skip = 2)
DEG_heat<- DEG_heat[,2]
names(DEG_heat)<-"Gene_ID"

DEG_heat<-DEG_heat%>%
  mutate(Gene_ID= gsub("CA.PGAv.1.6.","",Gene_ID))
  

top30_DEG_heat<- as.matrix(DEG_heat[1:30,])

# subsetting the countData according to DEG of heat stress
#top_heat_countData<- heat_logcount[top30_DEG_heat,]

# function to calculate Z-score
cal_z_score<- function(x){(x-mean(x))/sd(x)}

# calculate Zscore of top_30 heat countData

heat_Zscore<- t(apply(heat_logcount, 1 ,cal_z_score))
heat_Zscore<- heat_Zscore[top30_DEG_heat,]

# preparation for heatmap
# annotations
heat_annotation<- data.frame(Sample=colnames(heat_logcount))%>%
 mutate(Group= case_when(
    str_detect(Sample,"mock")~"mock",
    str_detect(Sample, "heat")~"heat"
  ))%>%
 column_to_rownames(var = "Sample")

annotation_color<- list(Group= c('mock' = 'gray47',
                                   'heat' = 'red'))

p1<- pheatmap(heat_Zscore,
         name= "Z-score",
         legend = TRUE,
         annotation_col = heat_annotation,
         annotation_colors = annotation_color,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "Heat",
         border_color = NA)


#----------------------------------------------------------------------#

####------------------------ Cold ---------------------------------#####

#----------------------------------------------------------------------#

# Steps for normalization
# step1: form DGEList
y_cold<- DGEList(cold_df)
y_cold$samples

#step2:  forming group 
cold_colData$group <- with(cold_colData, paste(stress,time,sep="."))
cold_colData$group<- factor(cold_colData$group)
cold_colData

# step3: filter the lowly expressed genes
cold_keep<- filterByExpr(y_cold, group = cold_colData$group)
y_cold<- y_cold[cold_keep, ,keep.lib.size= FALSE]
dim(y_cold)

# step4: normalization by TMM method
cold_TMM<- calcNormFactors(y_cold, 
                           group= cold_colData$group,
                           method = "TMM")
cold_logcount<- cpm(cold_TMM, log = TRUE)


# load the DEG data
# DEG of cold stress
library(readxl)
DEG_cold<- read_excel("DEG result for each stress/DEG result with cold stress.xlsx", skip = 2)
DEG_cold<- DEG_cold[,2]
names(DEG_cold)<-"Gene_ID"

DEG_cold<-DEG_cold%>%
  mutate(Gene_ID= gsub("CA.PGAv.1.6.","",Gene_ID))


top30_DEG_cold<- as.matrix(DEG_cold[1:30,])

# subsetting the countData according to DEG of cold stress
#top_cold_countData<- cold_logcount[top30_DEG_cold,]

# function to calculate Z-score
cal_z_score<- function(x){(x-mean(x))/sd(x)}

# calculate Zscore of top_30 cold countData

cold_Zscore<- t(apply(cold_logcount, 1 ,cal_z_score))
cold_Zscore<- cold_Zscore[top30_DEG_cold,]

# preparation for heatmap
# annotations
cold_annotation<- data.frame(Sample=colnames(cold_logcount))%>%
  mutate(Group= case_when(
    str_detect(Sample,"mock")~"mock",
    str_detect(Sample, "cold")~"cold"
  ))%>%
  column_to_rownames(var = "Sample")
cold_annotation

cold_annotation_color<- list(Group= c('mock' = 'gray47', 'cold' = 'royalblue'))

p2<- pheatmap(cold_Zscore,
              legend = TRUE,
              legend_labels = c("Z-score"),
              annotation_col = cold_annotation,
              annotation_colors = cold_annotation_color,
              cluster_rows = TRUE,
              cluster_cols = FALSE,
              main = "cold",
              border_color = NA)


#----------------------------------------------------------------------#

####------------------------ Mannitol -----------------------------#####

#----------------------------------------------------------------------#

# Steps for normalization
# step1: form DGEList
y_man<- DGEList(man_df)
y_man$samples

#step2:  forming group 
man_colData$group <- with(man_colData, paste(stress,time,sep="."))
man_colData$group<- factor(man_colData$group)
man_colData

# step3: filter the lowly expressed genes
man_keep<- filterByExpr(y_man, group = man_colData$group)
y_man<- y_man[man_keep, ,keep.lib.size= FALSE]
dim(y_man)

# step4: normalization by TMM method
man_TMM<- calcNormFactors(y_man, 
                           group= man_colData$group,
                           method = "TMM")
man_logcount<- cpm(man_TMM, log = TRUE)


# load the DEG data
# DEG of man stress
library(readxl)
DEG_man<- read_excel("DEG result for each stress/DEG result with man stress.xlsx", skip = 2)
DEG_man<- DEG_man[,2]
names(DEG_man)<-"Gene_ID"

DEG_man<-DEG_man%>%
  mutate(Gene_ID= gsub("CA.PGAv.1.6.","",Gene_ID))


top30_DEG_man<- as.matrix(DEG_man[1:30,])

# subsetting the countData according to DEG of man stress
#top_man_countData<- man_logcount[top30_DEG_man,]

# function to calculate Z-score
cal_z_score<- function(x){(x-mean(x))/sd(x)}

# calculate Zscore of top_30 man countData

man_Zscore<- t(apply(man_logcount, 1 ,cal_z_score))
man_Zscore<- man_Zscore[top30_DEG_man,]


# preparation for heatmap
# annotations
man_annotation<- data.frame(Sample=colnames(man_logcount))%>%
  mutate(Group= case_when(
    str_detect(Sample,"mock")~"mock",
    str_detect(Sample, "man")~"man"
  ))%>%
  column_to_rownames(var = "Sample")
man_annotation

man_annotation_color<- list(Group= c('mock' = 'gray47', 'man' = '#feb81a'))


p3<- pheatmap(man_Zscore,
         legend = TRUE,
         annotation_col = man_annotation,
         annotation_colors = man_annotation_color,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "Man",
         border_color = NA)


#----------------------------------------------------------------------#

####-------------------------- NaCl -------------------------------#####

#----------------------------------------------------------------------#


# Steps for normalization
# step1: form DGEList
y_NaCl<- DGEList(NaCl_df)
y_NaCl$samples

#step2:  forming group 
NaCl_colData$group <- with(NaCl_colData, paste(stress,time,sep="."))
NaCl_colData$group<- factor(NaCl_colData$group)
NaCl_colData

# step3: filter the lowly expressed genes
NaCl_keep<- filterByExpr(y_NaCl, group = NaCl_colData$group)
y_NaCl<- y_NaCl[NaCl_keep, ,keep.lib.size= FALSE]
dim(y_NaCl)

# step4: normalization by TMM method
NaCl_TMM<- calcNormFactors(y_NaCl, 
                          group= NaCl_colData$group,
                          method = "TMM")
NaCl_logcount<- cpm(NaCl_TMM, log = TRUE)


# load the DEG data
# DEG of NaCl stress
library(readxl)
DEG_NaCl<- read_excel("DEG result for each stress/DEG result with NaCl stress.xlsx", skip = 2)
DEG_NaCl<- DEG_NaCl[,2]
names(DEG_NaCl)<-"Gene_ID"

DEG_NaCl<-DEG_NaCl%>%
  mutate(Gene_ID= gsub("CA.PGAv.1.6.","",Gene_ID))


top30_DEG_NaCl<- as.matrix(DEG_NaCl[1:30,])

# subsetting the countData according to DEG of NaCl stress
#top_NaCl_countData<- NaCl_logcount[top30_DEG_NaCl,]

# function to calculate Z-score
cal_z_score<- function(x){(x-mean(x))/sd(x)}

# calculate Zscore of top_30 NaCl countData

NaCl_Zscore<- t(apply(NaCl_logcount, 1 ,cal_z_score))
NaCl_Zscore<- NaCl_Zscore[top30_DEG_NaCl,]

# preparation for heatmap
# annotations
NaCl_annotation<- data.frame(Sample=colnames(NaCl_logcount))%>%
  mutate(Group= case_when(
    str_detect(Sample,"mock")~"mock",
    str_detect(Sample, "NaCl")~"NaCl"
  ))%>%
  column_to_rownames(var = "Sample")
NaCl_annotation

NaCl_annotation_color<- list(Group= c('mock' = 'gray47', 'NaCl' = '#65a947'))

p4<- pheatmap(NaCl_Zscore,
         legend = TRUE,
         legend_labels = c("Z-score"),
         annotation_col = NaCl_annotation,
         annotation_colors = NaCl_annotation_color,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "NaCl",
         border_color = NA)


grid.arrange(p1$gtable, p2$gtable, 
             p3$gtable, p4$gtable, 
             nrow=2, ncol=2)





