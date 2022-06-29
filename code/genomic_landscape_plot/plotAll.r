
args<-commandArgs(TRUE)

## I/O
lineWidth <- 2
baseSize <- 45

##==================================================================
## Load libraries
##==================================================================
write("--------Loading libraries",stdout())
library(ggplot2)
library(zoo)
library(gtable)
library(gridExtra)
library(ggpubr)
library(gridGraphics)
library(cowplot)

##==================================================================
## Chromosome ordering
##==================================================================
chromosomeOrder <- c("chr1H", "chr2H", "chr3H", "chr4H", "chr5H", "chr6H", "chr7H", "chrUn")

##==================================================================
## Coordinates
##==================================================================
write("--------Read and plot coordinates", stdout())
coor <- read.table(args[1], header=T, sep="\t")
coor$chr <- factor(coor$CHROM,levels=chromosomeOrder)

xlab <- c(0,200,400,600)

coorPlot <- ggplot(coor) +
geom_segment(aes(x=coor$START/1000, y=1, xend=coor$END/1000, yend=1),size=12,color=coor$COLOR) +
scale_x_continuous(position = "top",labels=xlab,breaks=c(0,200000,400000,600000)) +
facet_wrap(~ chr,nrow=1) +
theme_bw(base_size = baseSize) +
theme(strip.background=element_rect(fill="white"), panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank()) +
theme(axis.title=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()) +
theme(strip.text.x=element_blank())+
theme(legend.position = "none")

##==================================================================
## Pi plot
##==================================================================
write("--------Read and plot pi file", stdout())
data <- read.table(args[2], header=T, sep="\t")
data$chr <- factor(data$CHROM,levels=chromosomeOrder)
dataPiSubset <- data[grep ("chrUn", data$chr, invert=TRUE), ]

rollAvgWindowSize <- as.numeric(10000)

piPlot<- ggplot(dataPiSubset, aes(POS/1000, PI)) + 
geom_line(size = lineWidth, colour="#b2182b",aes(x=dataPiSubset$POS/1000, y=rollmean(dataPiSubset$spontaneums, rollAvgWindowSize, na.pad=TRUE))) + 
geom_line(size = lineWidth, colour="#ef8a62",aes(x=dataPiSubset$POS/1000, y=rollmean(dataPiSubset$landraces, rollAvgWindowSize, na.pad=TRUE))) + 
geom_line(size = lineWidth, colour="#2166ac",aes(x=dataPiSubset$POS/1000, y=rollmean(dataPiSubset$cultivars, rollAvgWindowSize, na.pad=TRUE))) +
facet_wrap(~ dataPiSubset$chr,nrow=1) + 
scale_y_continuous(breaks=c(0,0.05,0.1),name="a")+
theme_classic(base_size = baseSize) +
theme(strip.background=element_rect(fill="white"), panel.border=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank()) +
theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank())+
theme(axis.title.y=element_text(angle = 0,face="bold"))+
theme(strip.text.x=element_blank())

##==================================================================
## Fst plot
##==================================================================

## Cult vs Spot
write("--------Read and plot Fst file", stdout())
data <- read.table(args[3], header=T, sep="\t")
data$chr <- factor(data$CHROM,levels=chromosomeOrder)
dataFst1Subset <- data[grep ("chrUn", data$chr, invert=TRUE), ]

fstPlot1<- ggplot(dataFst1Subset, aes(POS/1000, WEIR_AND_COCKERHAM_FST)) +
geom_point(size=1,color=ifelse(dataFst1Subset$WEIR_AND_COCKERHAM_FST>=0.8,"red","black")) +
facet_wrap(~ chr,nrow=1) +
scale_y_continuous(breaks=c(0,0.5,1),name="b",limits=c(0,1))+
theme_classic(base_size = baseSize) +
theme(strip.background=element_blank(), panel.border=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank()) +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
theme(axis.title.y=element_text(angle = 0,face="bold"))+
theme(strip.text.x=element_blank())

## Land vs Spot
write("--------Read and plot Fst file", stdout())
data <- read.table(args[4], header=T, sep="\t")
data$chr <- factor(data$CHROM,levels=chromosomeOrder)
dataFst2Subset <- data[grep ("chrUn", data$chr, invert=TRUE), ]

fstPlot2<- ggplot(dataFst2Subset, aes(POS/1000, WEIR_AND_COCKERHAM_FST)) +
geom_point(size=1,color=ifelse(dataFst2Subset$WEIR_AND_COCKERHAM_FST>=0.8,"red","black")) +
facet_wrap(~ chr,nrow=1) +
scale_y_continuous(breaks=c(0,0.5,1),name="c",limits=c(0,1))+
theme_classic(base_size = baseSize) +
theme(strip.background=element_blank(), panel.border=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank()) +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
theme(axis.title.y=element_text(angle = 0,face="bold"))+
theme(strip.text.x=element_blank())

## Cult vs Land
write("--------Read and plot Fst file", stdout())
data <- read.table(args[5], header=T, sep="\t")
data$chr <- factor(data$CHROM,levels=chromosomeOrder)
dataFst3Subset <- data[grep ("chrUn", data$chr, invert=TRUE), ]

fstPlot3<- ggplot(dataFst3Subset, aes(POS/1000, WEIR_AND_COCKERHAM_FST)) +
geom_point(size=1,color=ifelse(dataFst3Subset$WEIR_AND_COCKERHAM_FST>=0.8,"red","black")) +
facet_wrap(~ chr,nrow=1) +
scale_y_continuous(breaks=c(0,0.5,1),name="d",limits=c(0,1))+
theme_classic(base_size = baseSize) +
theme(strip.background=element_blank(), panel.border=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank()) +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
theme(axis.title.y=element_text(angle = 0,face="bold"))+
theme(strip.text.x=element_blank())


##==================================================================
## Haplotype block
##==================================================================
write("--------Read and plot haplotype block", stdout())
blocks <- read.table(args[6],sep="\t",header=T)
blocks$CHR <- factor(blocks$CHR,levels=chromosomeOrder)
blocks <- blocks[grep ("chrUn", blocks$CHR, invert=TRUE), ]

blockPlot <- ggplot(blocks, colour=CHR) + 
geom_segment(aes(x=blocks$BP1, y=blocks$GROUP, xend=blocks$BP2, yend=blocks$GROUP),size=15) + 
facet_wrap(~ CHR,nrow=1) +
ylab("e")+
theme_bw(base_size = baseSize) +
theme(strip.background=element_blank(), panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank()) +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
theme(axis.title.y=element_text(angle = 0,face="bold"))+
theme(strip.text.x=element_blank())

##==================================================================
## Tajima's D
##==================================================================
write("--------Read and plot Tajima's D", stdout())
data <- read.table(args[7], header=T, na.string="nan")
data$chr <- factor(data$CHROM,levels=chromosomeOrder)
dataTajSubset <- data[grep ("chrUn", data$chr, invert=TRUE), ]

rollAvgWindowSize <- as.integer(1000)

tajPlot<- ggplot(dataTajSubset, aes(BIN_START/1000, PI)) + 
geom_point(size = 0.3, colour="#d18884", aes(x=dataTajSubset$BIN_START/1000, y=dataTajSubset$spontaneums), na.rm=TRUE, alpha=0.8) +
geom_point(size = 0.3, colour="#f0b299", aes(x=dataTajSubset$BIN_START/1000, y=dataTajSubset$landraces),   na.rm=TRUE, alpha=0.8) +
geom_point(size = 0.3, colour="#7996b3", aes(x=dataTajSubset$BIN_START/1000, y=dataTajSubset$cultivars),   na.rm=TRUE, alpha=0.8) +
geom_line(size = lineWidth+1, colour="#b2182b", aes(x=dataTajSubset$BIN_START/1000, y=rollapply(as.numeric(as.character(dataTajSubset$spontaneums)), rollAvgWindowSize, mean, na.rm=TRUE, partial=TRUE, fill = NA)) ) +
geom_line(size = lineWidth+1, colour="#ef8a62", aes(x=dataTajSubset$BIN_START/1000, y=rollapply(as.numeric(as.character(dataTajSubset$landraces)), rollAvgWindowSize, mean, na.rm=TRUE, partial=TRUE, fill = NA)) ) + 
geom_line(size = lineWidth+1, colour="#2166ac", aes(x=dataTajSubset$BIN_START/1000, y=rollapply(as.numeric(as.character(dataTajSubset$cultivars)), rollAvgWindowSize, mean, na.rm=TRUE, partial=TRUE, fill = NA)) ) + 
#geom_smooth(size=lineWidth, colour="#b2182b", aes(x=dataTajSubset$BIN_START/1000, y=dataTajSubset$spontaneums))+
#geom_smooth(size=lineWidth, colour="#ef8a62", aes(x=dataTajSubset$BIN_START/1000, y=dataTajSubset$landraces))+
#geom_smooth(size=lineWidth, colour="#2166ac", aes(x=dataTajSubset$BIN_START/1000, y=dataTajSubset$cultivars))+
geom_hline(aes(yintercept=1.19227))+
geom_hline(aes(yintercept=-2.04366))+
facet_wrap(~ dataTajSubset$chr,nrow=1) + 
ylab("f")+
theme_gray(base_size = baseSize) +
theme(strip.background=element_blank(), panel.border=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank()) +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
theme(axis.title.y=element_text(angle = 0,face="bold"))+
theme(strip.text.x=element_blank())

##==================================================================
## Mu (selective sweeps)
##==================================================================
write("--------Read and plot RAiSD Mu", stdout())
data <- read.table(args[8], header=T, na.string="nan")
data$chr <- factor(data$CHR,levels=chromosomeOrder)
dataMuSubset <- data[grep ("chrUn", data$chr, invert=TRUE), ]

## Cultivar plot
dataCultSubset <- dataMuSubset[grep("cultivars",dataMuSubset$GROUP),]

muPlot3<- ggplot(dataCultSubset, aes(x=LOC/1000, y=MU)) +
geom_point(size = lineWidth+1, color=ifelse(dataCultSubset$MU>=4.566e-05,"red","black")) +
facet_wrap(~ dataCultSubset$chr,nrow=1) +
ylab("i")+
theme_gray(base_size = baseSize) +
theme(strip.background=element_blank(), panel.border=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank()) +
theme(axis.title.x=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) +
theme(strip.text.x=element_blank())+
theme(axis.title.y=element_text(angle = 0,face="bold"))+
theme(legend.position = "none")

## Landrace plot
dataLandSubset <- dataMuSubset[grep("landraces",dataMuSubset$GROUP),]

muPlot2<- ggplot(dataLandSubset, aes(x=LOC/1000, y=MU)) +
geom_point(size = lineWidth+1, color=ifelse(dataLandSubset$MU>=1.993e-05,"red","black")) +
facet_wrap(~ dataLandSubset$chr,nrow=1) +
ylab("h")+
theme_gray(base_size = baseSize) +
theme(strip.background=element_blank(), panel.border=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank()) +
theme(axis.title.x=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) +
theme(strip.text.x=element_blank())+
theme(axis.title.y=element_text(angle = 0,face="bold"))+
theme(legend.position = "none")

## Spontaneums plot
dataSponSubset <- dataMuSubset[grep("spontaneums",dataMuSubset$GROUP),]

muPlot1<- ggplot(dataSponSubset, aes(x=LOC/1000, y=MU)) +
geom_point(size = lineWidth+1, color=ifelse(dataSponSubset$MU>=1.267e-06,"red","black")) +
facet_wrap(~ dataSponSubset$chr,nrow=1) +
ylab("g")+
theme_gray(base_size = baseSize) +
theme(strip.background=element_blank(), panel.border=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank()) +
theme(axis.title.x=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) +
theme(strip.text.x=element_blank())+
theme(axis.title.y=element_text(angle = 0,face="bold"))+
theme(legend.position = "none")

##==================================================================
## Domestication genes
##==================================================================
write("--------Read and plot domestication genes", stdout())
dom <- read.table(args[9], header=T, sep="\t")
dom$chr <- factor(dom$Chrom,levels=chromosomeOrder)

domPlot <- ggplot(dom) +
geom_segment(aes(x=dom$Start/1000, y=1, xend=dom$End/1000, yend=1),size=20,color=dom$Colour) +
facet_wrap(~ chr,nrow=1) +
ylab("j")+
theme_bw(base_size = baseSize) +
theme(strip.background=element_rect(fill="white"), panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank())+
theme(axis.title.x=element_blank(), axis.ticks=element_blank(), axis.text=element_blank()) +
theme(strip.text.x=element_blank())+
theme(axis.title.y=element_text(angle = 0,face="bold"))+
theme(legend.position = "none")

##==================================================================
## Plotting the panel
##==================================================================
write("--------Running ggplotGrob coordinates",stdout())
g0 <- ggplotGrob(coorPlot)
write("--------Running ggplotGrob pi plot",stdout())
g1 <- ggplotGrob(piPlot)
write("--------Running ggplotGrob fst plot 1",stdout())
g2.1 <- ggplotGrob(fstPlot1)
write("--------Running ggplotGrob fst plot 2",stdout())
g2.2 <- ggplotGrob(fstPlot2)
write("--------Running ggplotGrob fst plot 3",stdout())
g2.3 <- ggplotGrob(fstPlot3)
write("--------Running ggplotGrob block plot",stdout())
#g3 <- ggplotGrob(blockPlot)
#write("--------Running ggplotGrob tajima's D plot",stdout())
#g4 <- ggplotGrob(tajPlot)
#write("--------Running ggplotGrob RAiSD Mu plot1",stdout())
#g5.1 <- ggplotGrob(muPlot1)
#write("--------Running ggplotGrob RAiSD Mu plot2",stdout())
#g5.2 <- ggplotGrob(muPlot2)
#write("--------Running ggplotGrob RAiSD Mu plot3",stdout())
#g5.3 <- ggplotGrob(muPlot3)
#write("--------Running ggplotGrob domestication genes",stdout())
#g6 <- ggplotGrob(domPlot)


write("--------Creating panel",stdout())
## cowplot method. easier to specify row height
#plotPanel <- plot_grid(g0,g1,g2.1,g2.2,g2.3,g3,g4,g5.1,g5.2,g5.3,g6,align="v",nrow=11,rel_heights=c(0.5,0.8,1,1,1,0.65,0.8,0.5,0.5,0.5,0.5))
plotPanel <- plot_grid(g0,g1,g2.1,g2.2,g2.3,align="v",nrow=5,rel_heights=c(0.8,1,1,1,1))

write("--------Export to a PNG file", stdout())
## Original figure
#png(args[10],2000,1300)
## New figure
png(args[10],3000,1200)
grid.draw(plotPanel)
dev.off()
write("--------Done",stdout())
