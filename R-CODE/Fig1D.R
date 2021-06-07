
## Author: Mostafa Rahnama 
## Date: 1/1/2021
## purpose: plotting TAS blast results and random sampling data (Fig1D) 

library(readxl)
library(ggplot2)
library(tidyverse)
library(ggh4x)
####### 
strains <- read.csv("~/TAS_BLAST_Mastersheet_Final.csv")

strains_l <- gather(strains, subject_isolate, hits, Bm8309_B1:Pg1054_St)
strains_l <- separate(data = strains_l, col = subject_isolate, into = c("subject_isolate", "subject_clade"), sep = "_")   ###### be careful not to have '_' in strain name
strains_l$strain <- sub("(.*)_.*", "\\1",  strains_l$consensus_ID) 

# remove columns with zero hits
strains_l <- subset(strains_l, hits != 0 )
strains_l2 <- ddply(strains_l, c("strain","query_clade","consensus_ID"), summarise, Clade_hits= length(unique(subject_clade)), 
                                                                    strain_hits= length(unique(subject_isolate)), aver_hits= mean(hits))
strains_l2$color <- ifelse(as.numeric(strains_l2$Clade_hits) == 1, 'Lineage specific', 'Lineage non-specific')
strains_l2$levels <- paste("telConsensus")

########## Reading random sampling data
Random <- read.csv("~/Random_Sampling280bp_Final.csv")
Random2 <- ddply(Random, c("strain", "query_clade","consensus_ID"), summarise, Clade_hits= length(unique(subject_clade)), 
                   strain_hits= length(unique(subject_isolate)), aver_hits= mean(hits))
Random2$color <- ifelse(as.numeric(Random2$Clade_hits) == 1, 'Lineage specific (Random)', 'Lineage non-specific (Random)')
Random2$levels <- paste("R")   #Random
### mix both table
Fianl_df <- rbind(Random2, strains_l2)

Fianl_df$query_clade <- ifelse(Fianl_df$query_clade == "B1", "U1", as.character(Fianl_df$query_clade))
Fianl_df$query_clade <- ifelse(Fianl_df$query_clade == "B2", "U2", as.character(Fianl_df$query_clade))
Fianl_df$query_clade <- ifelse(Fianl_df$query_clade == "B3", "U3", as.character(Fianl_df$query_clade))
Fianl_df$query_clade <- ifelse(Fianl_df$query_clade == "B4", "U4", as.character(Fianl_df$query_clade))

Fianl_df$query_clade <- factor(Fianl_df$query_clade, 
                               levels = c("C1","C2","Ce","E1","E2","E3","Ec","Er","H","L1","L2",
                                          "L3","Le","Lu","O","P","S","St","T","Pg","Pp","Pu","U1","U2","U3","U4"))

Fianl_df$levels <- ifelse(Fianl_df$levels == "telConsensus", "TAS", (Fianl_df$levels))
Fianl_df$class <- factor(Fianl_df$levels, levels = c("TAS", "R"))

# function for computing mean, DS, max and min values
min.mean.sd.max <- function(x) {
  r <- c(min(x), quantile(x,0.25), median(x), quantile(x,0.75), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r}


  ggplot(Fianl_df, aes(x = query_clade, y = as.numeric(strain_hits)))+
    stat_summary(fun.data = min.mean.sd.max, geom = "boxplot") + 
    geom_jitter(aes(color = color),position=position_jitter(width=.2), size=3) +   
    labs(title = "", x = "\nPopulation", y= "Number of strains containing query sequence\n", color = "")+
    scale_color_manual(values = c("Lineage specific" = "blue","Lineage non-specific" = "red",
                                  "Lineage specific (Random)" = "chartreuse4","Lineage non-specific (Random)" = "#ff8303"))+
    facet_nested(. ~ query_clade+class, scales = "free", space = "free", switch ="x",nest_line = FALSE, bleed = FALSE)+   #cols = vars(Chr) 
    theme(#plot.title = element_text(hjust = 1.5),
      legend.title=element_text(face="bold", size=11), 
      legend.text=element_text(face="bold", size=19),
      legend.direction = "horizontal",
      legend.key = element_rect(fill = "white"),
      legend.key.size = unit(1, "cm"),
      legend.position = "bottom", 
      legend.justification="bottom",
      axis.title.x = element_text(face="bold", size = 20),  # , margin = margin(1, 0, 40, 0), vjust = -3
      axis.title.y = element_text(face="bold", size = 20),  
      axis.text.y=element_text(size=14, face="bold", colour = "black", angle=0),
      axis.text.x = element_blank(),
      axis.line.y = element_line(colour = "black", size = 0.2),
      strip.text.x = element_text(size=14, face = "bold", colour = "black", angle=90), #, margin = margin(0.15,0,0.15,0, "cm")
      strip.background=element_rect(color="grey30", fill="grey90"),
      panel.spacing=unit(0.5,"lines"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(), 
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.background = element_rect(fill = "white", colour = "white", size = 0.4),
      panel.border = element_rect(colour ="white", fill=NA, size=0.1),
      axis.ticks.length.y = unit(0.1,"cm"),
      axis.ticks.y = element_line(size = 0.5),
      axis.ticks.x = element_blank()
    )
  + guides(color = guide_legend(override.aes = list(size = 5)))



