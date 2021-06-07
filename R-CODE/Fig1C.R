
## Authour: Mostafa Rahnama 
## Date: 1/1/2021
## purpose: plotting TJ data (Fig1C) 

library(readxl)
library(ggplot2)
library(tidyverse)
library(dplyr)

strains4_l <- read.csv("~/Google Drive/1.Forntier_Paper/1.Frontier_manusript/FINAL_VERSION/Code/TJ_Final.csv")
strains5_l <- ddply(strains4_l, c("query_clade","consensus_ID", "strain"), summarise, Clade_hits= length(unique(subject_clade)), 
                   strain_hits= length(unique(subject_isolate)), aver_hits= mean(hits))

# function for computing mean, DS, max and min values
min.mean.sd.max <- function(x) {
  r <- c(min(x), quantile(x,0.25), median(x), quantile(x,0.75), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r}

strains_Sum <- ddply(strains5_l, c("query_clade"), summarise, mean= mean(strain_hits), Medain= median(strain_hits))

strains5_l$color <- ifelse(as.numeric(strains5_l$Clade_hits) == 1, 'Lineage specific', 'Lineage non-specific')
strains5_l$query_clade <- ifelse(strains5_l$query_clade == "B1", "U1", as.character(strains5_l$query_clade))
strains5_l$query_clade <- ifelse(strains5_l$query_clade == "B2", "U2", as.character(strains5_l$query_clade))
strains5_l$query_clade <- ifelse(strains5_l$query_clade == "B3", "U3", as.character(strains5_l$query_clade))
strains5_l$query_clade <- ifelse(strains5_l$query_clade == "B4", "U4", as.character(strains5_l$query_clade))

strains5_l$query_clade <- factor(strains5_l$query_clade, 
                                levels = c("C1","C2","Ce","E1","E2","E3","Ec","Er","H","L1","L2",
                                           "L3","Le","Lu","O","P","S","St","T","U1","U2","U3","U4","Pg","Pp","Pu"))
strains5_l <- strains5_l[complete.cases(strains5_l),]

ggplot(strains5_l, aes(x = query_clade, y = as.numeric(strain_hits)))+
    stat_summary(fun.data = min.mean.sd.max, geom = "boxplot") + 
    geom_jitter(aes(color = color),position=position_jitter(width=.3, h=0.1), size=3) +   
    labs(title = "", x = "\nPopulation", y= "Number of strains in each consensus \n", color = "")+
    scale_color_manual(values = c("Lineage specific" = "blue","Lineage non-specific" = "red"))+
    coord_cartesian(ylim = c(0,10)) + 
    scale_y_continuous(breaks = c(0,1, 2,3,4,5,6,7,8,9,10)) +
    theme(#plot.title = element_text(hjust = 1.5),

      axis.title.x = element_text(face="bold", size = 20),  # , margin = margin(1, 0, 40, 0), vjust = -3
      axis.title.y = element_text(face="bold", size = 20),  
      axis.text.y=element_text(size=14, face="bold", colour = "black", angle=0, vjust=1, hjust = 0.5),
      axis.text.x = element_text(size=14, face="bold", colour = "black", angle=0, vjust=1, hjust = 0.5),
      strip.text.x = element_text(size=14, face = "bold", colour = "black", angle=90), #, margin = margin(0.15,0,0.15,0, "cm")
      strip.background=element_rect(color="grey30", fill="grey90"),
      panel.spacing=unit(0.5,"lines"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(), 
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.background = element_rect(fill = "white", colour = "black", size = 0.2),
      panel.border = element_rect(colour ="white", fill=NA, size=0.1),
      #ggh4x.facet.nestline = element_line(linetype = 10),
      axis.ticks = element_blank()    #element_line(size=1),
      )



















