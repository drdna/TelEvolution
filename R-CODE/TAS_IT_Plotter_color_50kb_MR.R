library(lattice)
library(ggplot2)
library(dplyr)
#library(Cairo)
#library(gridExtra)
#library(grid)
#library(png)

# Specify chromosome lengths depending on reference used

#U233lookUp <- data.frame(chrLengths = c(5.500179, 7.853447, 7.955039, 7.641679, 5.440338, 3.716274, 1.631780), subject = c('Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'Chr6', 'Chr7'))
CD156lookUp <- data.frame(chrLengths = c(6.012563, 8.265001, 7.555746, 5.498619, 4.609210, 6.066300, 3.913369), subject = c('Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'Chr6', 'Chr7'))
#ArcadialookUp <- data.frame(chrLengths = c(5.987451, 8.805343, 6.770928, 7.569919, 4.820434, 6.3375610, 3.501709), subject = c('Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'Chr6', 'Chr7'))

#read in blast report
blastReport <- read.table("Desktop/TrimmedConsensus.CD156.BLASTunique_Redundant_removed_MR2.txt", header = FALSE)

#read in internal telomere positions
TelConsensReport <- read.csv("All_internal_Tels_50bpGap_relative.csv", header = TRUE)

colnames(blastReport) = c("query", "subject", "qlen", "slen", "percent", "length", "mm", "gaps", "qstart", "qend", "sstart", "send", "evalue", "score", "Duplication")

blastReport %>% gsub('_.*','', query)

colnames(TelConsensReport) = c("query", "type", "sstart", "send", "subject", "Seq", "preSeq", "postSeq")

#TelConsensReport <- filter(TelConsensReport, type == "Dimer")

TelConsensReport <- filter(TelConsensReport, !grepl("Dimer", type))

blastReport$defl = rep(0.035, nrow(blastReport))

blastReport$sstart <- blastReport$sstart/1000000

blastReport$send <- blastReport$send/1000000


blastReport$dir <- ifelse(blastReport$sstart < blastReport$send, 'F', 'R')


TelConsensReport$defl = rep(-0.035, nrow(TelConsensReport))

TelConsensReport$sstart <- TelConsensReport$sstart/1000000

TelConsensReport$send <- TelConsensReport$send/1000000

TelConsensReport$dir <- ifelse(grepl("CCCTAA", TelConsensReport$Seq), 'F', 'R')

#subjF <- rev(factor(subject))

#queryF <- factor(query)

levels(subjF)

# filter to plot single chromosome

#blastReport <- filter(blastReport, subject == "Chr6")
#TelConsensReport <- filter(TelConsensReport, subject == "Chr6")
#CD156lookUp <- filter(CD156lookUp, subject == "Chr6")
blastReport$Duplication <- factor(blastReport$Duplication, 
                                 levels = c("Y","N"))

#pdf(file="/Users/mostafa/Google Drive/1.Forntier_Paper/1.Frontier_manusript/DATA/FiG_Fig2_newest.pdf", encoding="MacRoman",8, 8)  #
ppi <- 300
tiff(file=paste("TJ/Fig2A",".tiff", sep = ''), width = 9*ppi, height= 14*ppi, res = ppi)

print(
ggplot() + 
  geom_segment(data = CD156lookUp, aes(x=0, xend=chrLengths, y=-0.004, yend=-0.004, size = 1, colour= "segment")) +
  geom_jitter(data = blastReport, aes(x=sstart, y = defl, shape = factor(dir), colour = Duplication), size = 4, height=0.025) + 
  scale_shape_manual(values=c("\u25BA", "\u25C4")) +
  scale_colour_manual(values = c("Y"="Red", "N"=" Black", "Dimer"="green", "Trimer"="orange", "Tetramer"="blue","segment"="black")) +
  geom_jitter(data = TelConsensReport, aes(x=sstart, y = defl, shape = factor(dir), colour = type), size = 4, height=0.025) + 
  scale_shape_manual(values=c("\u25BA", "\u25C4")) +
 
 #ylim(-0.05, 0.06) + # xlim(6.01, 6.06) +

  facet_wrap(~ subject, ncol = 1) +
 #  scale_shape_identity() +
  theme(strip.text.x = element_text(size = 12, face = "bold")) + # x-out for un-annotated plot
  #theme(strip.text.x = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank()) +
  xlab("Chromosome position (kb)") +   # x-out for un-annotated plot
  theme(legend.position = "none") + # x-out for un-annotated plot
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) #+
  #ggtitle("TAS positions in CD156 genome") # x-out for un-annotated plot
)
dev.off()



### subset of chr 6
CD156lookUp_sub <- data.frame(chrLengths = c(6.06), subject = c( 'Chr6'))
blastReport_sub <- subset(blastReport, subject == "Chr6")  %>% subset(sstart > 6)
TelConsensReport_sub <- subset(TelConsensReport, subject == "Chr6")  %>% subset(sstart > 6)

#blastReport_sub$dir <- ifelse(blastReport_sub$dir == "R", 'F', 'R')
#TelConsensReport_sub$dir <- ifelse(TelConsensReport_sub$dir == "R", 'R', 'R')
ppi <- 300
#tiff(file=paste("/Users/mostafa/Google Drive/1.Forntier_Paper/1.Frontier_manusript/DATA/Fig2B",".tiff", sep = ''),
     width = 6*ppi, height= 3.2*ppi, res = ppi)

#pdf(file=paste("/Users/mostafa/Google Drive/1.Forntier_Paper/1.Frontier_manusript/DATA/Fig2B", 
#               ".pdf", sep = ''), 8, 8, useDingbats=FALSE)
print(
ggplot() + 
  geom_segment(data = CD156lookUp_sub, aes(x=6, xend=chrLengths, y=-0.004, yend=-0.004, size = 1, colour= "segment")) +
  geom_jitter(data = blastReport_sub, aes(x=sstart, y = defl, shape = factor(dir), colour = Duplication), size = 4, height=0.025) + 
  scale_shape_manual(values=c("\u25C4")) +
  scale_colour_manual(values = c("Y"="Red", "N"=" Black", "Dimer"="green", "Trimer"="orange", "Tetramer"="blue"
                                 ,"segment"="black"))+
  geom_jitter(data = TelConsensReport_sub, aes(x=sstart, y = defl, shape = factor(dir), colour = type), size = 4, height=0.025) + 
  scale_shape_manual(values=c("\u25C4")) +
  
  ylim(-0.06, 0.06) + 
  xlim(6, NA) +
  
  facet_wrap(~ subject, ncol = 1) +
  #  scale_shape_identity() +
  theme(strip.text.x = element_text(size = 12, face = "bold")) + # x-out for un-annotated plot
  #theme(strip.text.x = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank()) +
  xlab("Chromosome position (kb)") +   # x-out for un-annotated plot
  theme(legend.position = "none") + # x-out for un-annotated plot
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) #+
#ggtitle("TAS positions in CD156 genome") # x-out for un-annotated plot
)
#dev.off()


