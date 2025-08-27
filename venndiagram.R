
getwd()
setwd("C:/Users/ASUS/OneDrive/Desktop/Chili/DEG result for each stress")

# load the libraries
library(tidyverse)
library(venn)
library(VennDiagram)
library(readxl)

# read the excel files
list.files()

# DEG of Heat stress
DEG_heat<- read_excel("DEG result with heat stress.xlsx", skip = 2)
DEG_heat<- DEG_heat[,2]
names(DEG_heat)<- "Gene_ID"

DEG_heat<-as.matrix(DEG_heat)

# DEG of Cold Stress
DEG_cold<- read_excel("DEG result with cold stress.xlsx", skip = 2)
DEG_cold<- DEG_cold[,2]
names(DEG_cold)<- "Gene_ID"
DEG_cold<-as.matrix(DEG_cold)

# DEG of Mannitol Stress
DEG_man<- read_excel("DEG result with man stress.xlsx", skip = 2)
DEG_man<- DEG_man[,2]
names(DEG_man)<- "Gene_ID"

DEG_man<- as.matrix(DEG_man)


# DEG of NaCl stress
DEG_NaCl<- read_excel("DEG result with NaCl stress.xlsx", skip = 2)
DEG_NaCl<- DEG_NaCl[,2]
names(DEG_NaCl)<- "Gene_ID"
DEG_NaCl<- as.matrix(DEG_NaCl)


# venn diagram
# venn data
venn_data<- list(Heat=DEG_heat,
                 NaCl= DEG_NaCl,
                 Cold=DEG_cold,
                 Man=DEG_man
                 )

# create venn colors
vennColors<- c('red',"#65a947",'royalblue','#feb81a')


venn_config<- venn.diagram(x= venn_data,
                           category.names = names(venn_data),
                           fill= vennColors,
                           filename = NULL,
                           imagetype = "png",
                           col="transparent",
                           alpha=0.5,
                           cex= 1.4,
                           fontfamily="sans",
                           cat.cex=1.2,
                           cat.fontfamily="sans",
                           margin=0.1,
                           cex.label=2,
                           lwd=2)

# Display the Venn diagram
grid.newpage()
grid.draw(venn_config)




