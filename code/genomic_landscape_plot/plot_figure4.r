args<-commandArgs(TRUE)

##===============================================
## Load libraries
##===============================================
write("--------Loading libraries",stdout())
library(ggplot2)
library(zoo)
library(gridExtra)
library(gridGraphics)
library(cowplot)

##===============================================
## I/O
##===============================================
write("--------Reading in data",stdout())
coor     <- read.table(args[1], header=T, sep="\t")
haplo    <- read.table(args[2], header=T, sep="\t")
dataMu   <- read.table(args[3], header=T, na.string="nan")
delsnps  <- read.table(args[4], header=T, sep="\t")
genes    <- read.table(args[5], header=T, sep="\t")
taj      <- read.table(args[6], header=T, na.string="nan")
outName  <- args[7]
baseSize=45

lineWidth <- as.integer(5)
chromosomeOrder <- c("chr1H", "chr2H", "chr3H", "chr4H", "chr5H", "chr6H", "chr7H", "chrUn")

##===============================================
## Plot zone3 coordinates
##===============================================
write("--------Plot zone3 scale",stdout())
coor$chr <- factor(coor$CHROM,levels=chromosomeOrder)

## Convert zone3 coordinate so plot x axis can start with 0
coor$START[coor$CHROM == "chr1H"] <- coor$START[coor$CHROM == "chr1H"] - 151200000
coor$START[coor$CHROM == "chr2H"] <- coor$START[coor$CHROM == "chr2H"] - 256000000
coor$START[coor$CHROM == "chr3H"] <- coor$START[coor$CHROM == "chr3H"] - 224000000
coor$START[coor$CHROM == "chr4H"] <- coor$START[coor$CHROM == "chr4H"] - 163200000
coor$START[coor$CHROM == "chr5H"] <- coor$START[coor$CHROM == "chr5H"] - 168800000
coor$START[coor$CHROM == "chr6H"] <- coor$START[coor$CHROM == "chr6H"] - 199200000
coor$START[coor$CHROM == "chr7H"] <- coor$START[coor$CHROM == "chr7H"] - 254400000
coor$END[coor$CHROM == "chr1H"] <- coor$END[coor$CHROM == "chr1H"] - 151200000
coor$END[coor$CHROM == "chr2H"] <- coor$END[coor$CHROM == "chr2H"] - 256000000
coor$END[coor$CHROM == "chr3H"] <- coor$END[coor$CHROM == "chr3H"] - 224000000
coor$END[coor$CHROM == "chr4H"] <- coor$END[coor$CHROM == "chr4H"] - 163200000
coor$END[coor$CHROM == "chr5H"] <- coor$END[coor$CHROM == "chr5H"] - 168800000
coor$END[coor$CHROM == "chr6H"] <- coor$END[coor$CHROM == "chr6H"] - 199200000
coor$END[coor$CHROM == "chr7H"] <- coor$END[coor$CHROM == "chr7H"] - 254400000

## Plotting
coorPlot <- ggplot(coor) +
	geom_segment(aes(x=coor$START/1000000, y=1, xend=coor$END/1000000, yend=1),size=8,color=coor$COLOR) +
	scale_x_continuous(position = "top",limits=c(0,243.2),breaks=c(0,50,100,150,200)) +
	facet_wrap(~ chr,nrow=1) +
	theme_bw(base_size = baseSize-5) +
	theme(strip.background=element_rect(fill="white"), panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank()) +
	theme(axis.title=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()) +
	theme(strip.text.x=element_blank())+
	theme(legend.position = "none")

##===============================================
## Plot number of gene haplotype
##===============================================
write("--------Plot gene haplotype count",stdout())

haplo$CHROM <- factor(haplo$CHROM,levels=chromosomeOrder)

## Convert zone3 coordinate
haplo$POS[haplo$CHROM == "chr1H"] <- haplo$POS[haplo$CHROM == "chr1H"] - 151200000
haplo$POS[haplo$CHROM == "chr2H"] <- haplo$POS[haplo$CHROM == "chr2H"] - 256000000
haplo$POS[haplo$CHROM == "chr3H"] <- haplo$POS[haplo$CHROM == "chr3H"] - 224000000
haplo$POS[haplo$CHROM == "chr4H"] <- haplo$POS[haplo$CHROM == "chr4H"] - 163200000
haplo$POS[haplo$CHROM == "chr5H"] <- haplo$POS[haplo$CHROM == "chr5H"] - 168800000
haplo$POS[haplo$CHROM == "chr6H"] <- haplo$POS[haplo$CHROM == "chr6H"] - 199200000
haplo$POS[haplo$CHROM == "chr7H"] <- haplo$POS[haplo$CHROM == "chr7H"] - 254400000

## Plotting
haploPlot<- ggplot(haplo, aes(POS/1000000, PI)) + 
	xlim(0, 243.2) +
	geom_point(stroke=4, size = lineWidth, colour="#b2182b", aes(x=haplo$POS/1000000, y=haplo$spontaneums), shape=0) +
	geom_point(stroke=4, size = lineWidth, colour="#ef8a62", aes(x=haplo$POS/1000000, y=haplo$landraces),   shape=1) + 
	geom_point(stroke=4, size = lineWidth, colour="#2166ac", aes(x=haplo$POS/1000000, y=haplo$cultivars),   shape=2) + 
	scale_y_continuous(breaks=c(0,2,4,6),limits=c(0,7))+
	facet_wrap(~ haplo$CHROM,nrow=1) +
	ylab("a")+
	theme_bw(base_size = baseSize) +
	theme(strip.background=element_rect(fill="white"), panel.border=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank()) +
	theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank())+
	theme(axis.title.y=element_text(angle = 0,face="bold",size=14))+
	theme(strip.text.x=element_blank()) 

##===============================================
## Plot selective sweep
##===============================================
write("--------Plot selective  sweeps",stdout())
dataMu$chr <- factor(dataMu$CHR,levels=chromosomeOrder)
dataMuSubset <- dataMu[grep ("chrUn", dataMu$chr, invert=TRUE), ]

## Keeping only the zone3 data
dataMuSubset <- subset(dataMuSubset, (chr=="chr1H" & LOC<=303200000)|(chr=="chr2H" & LOC<=401600000)|(chr=="chr3H" & LOC<=380000000)|(chr=="chr4H" & LOC<=406400000)|(chr=="chr5H" & LOC<=301600000)|(chr=="chr6H" & LOC<=316000000)|(chr=="chr7H" & LOC<=384000000))

## Convert zone3 coordinate
dataMuSubset$LOC[dataMuSubset$chr == "chr1H"] <- dataMuSubset$LOC[dataMuSubset$chr == "chr1H"] - 151200000
dataMuSubset$LOC[dataMuSubset$chr == "chr2H"] <- dataMuSubset$LOC[dataMuSubset$chr == "chr2H"] - 256000000
dataMuSubset$LOC[dataMuSubset$chr == "chr3H"] <- dataMuSubset$LOC[dataMuSubset$chr == "chr3H"] - 224000000
dataMuSubset$LOC[dataMuSubset$chr == "chr4H"] <- dataMuSubset$LOC[dataMuSubset$chr == "chr4H"] - 163200000
dataMuSubset$LOC[dataMuSubset$chr == "chr5H"] <- dataMuSubset$LOC[dataMuSubset$chr == "chr5H"] - 168800000
dataMuSubset$LOC[dataMuSubset$chr == "chr6H"] <- dataMuSubset$LOC[dataMuSubset$chr == "chr6H"] - 199200000
dataMuSubset$LOC[dataMuSubset$chr == "chr7H"] <- dataMuSubset$LOC[dataMuSubset$chr == "chr7H"] - 254400000

## Plotting
## Cultivar plot
dataCultSubset <- dataMuSubset[grep("cultivars",dataMuSubset$GROUP),]
muPlot3<- ggplot(dataCultSubset, aes(x=LOC/1000000, y=MU)) +
        xlim(0, 243.2) +
	geom_point(size = lineWidth+1, color=ifelse(dataCultSubset$MU>=4.566e-05,"red","black")) +
	facet_wrap(~ dataCultSubset$chr,nrow=1) +
	ylab("4")+
	theme_gray(base_size = baseSize) +
	theme(strip.background=element_blank(), panel.border=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank()) +
	theme(axis.title.x=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) +
	theme(strip.text.x=element_blank())+
	theme(legend.position = "none")

## Landrace plot
dataLandSubset <- dataMuSubset[grep("landraces",dataMuSubset$GROUP),]
muPlot2<- ggplot(dataLandSubset, aes(x=LOC/1000000, y=MU)) +
        xlim(0, 243.2) +
	geom_point(size = lineWidth+1, color=ifelse(dataLandSubset$MU>=1.993e-05,"red","black")) +
	facet_wrap(~ dataLandSubset$chr,nrow=1) +
	ylab("3")+
	theme_gray(base_size = baseSize) +
	theme(strip.background=element_blank(), panel.border=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank()) +
	theme(axis.title.x=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) +
	theme(strip.text.x=element_blank())+
	theme(legend.position = "none")

## Spontaneums plot
dataSponSubset <- dataMuSubset[grep("spontaneums",dataMuSubset$GROUP),]
muPlot1<- ggplot(dataSponSubset, aes(x=LOC/1000000, y=MU)) +
        xlim(0, 243.2) +
	geom_point(size = lineWidth+1, color=ifelse(dataSponSubset$MU>=1.267e-06,"red","black")) +
	facet_wrap(~ dataSponSubset$chr,nrow=1) +
	ylab("2")+
	theme_gray(base_size = baseSize) +
	theme(strip.background=element_blank(), panel.border=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank()) +
	theme(axis.title.x=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) +
	theme(strip.text.x=element_blank())+
	theme(legend.position = "none")

##===============================================
## Plotting deleterious SNPs
##===============================================
write("--------Plotting deleterious SNPs",stdout())
delsnps$chr <- factor(delsnps$Chrom,levels=chromosomeOrder)

## Convert zone 3 coordinate
delsnps$Start[delsnps$chr == "chr1H"] <- delsnps$Start[delsnps$chr == "chr1H"] - 151200000
delsnps$Start[delsnps$chr == "chr2H"] <- delsnps$Start[delsnps$chr == "chr2H"] - 256000000
delsnps$Start[delsnps$chr == "chr3H"] <- delsnps$Start[delsnps$chr == "chr3H"] - 224000000
delsnps$Start[delsnps$chr == "chr4H"] <- delsnps$Start[delsnps$chr == "chr4H"] - 163200000
delsnps$Start[delsnps$chr == "chr5H"] <- delsnps$Start[delsnps$chr == "chr5H"] - 168800000
delsnps$Start[delsnps$chr == "chr6H"] <- delsnps$Start[delsnps$chr == "chr6H"] - 199200000
delsnps$Start[delsnps$chr == "chr7H"] <- delsnps$Start[delsnps$chr == "chr7H"] - 254400000
delsnps$End[delsnps$chr == "chr1H"] <- delsnps$End[delsnps$chr == "chr1H"] - 151200000
delsnps$End[delsnps$chr == "chr2H"] <- delsnps$End[delsnps$chr == "chr2H"] - 256000000
delsnps$End[delsnps$chr == "chr3H"] <- delsnps$End[delsnps$chr == "chr3H"] - 224000000
delsnps$End[delsnps$chr == "chr4H"] <- delsnps$End[delsnps$chr == "chr4H"] - 163200000
delsnps$End[delsnps$chr == "chr5H"] <- delsnps$End[delsnps$chr == "chr5H"] - 168800000
delsnps$End[delsnps$chr == "chr6H"] <- delsnps$End[delsnps$chr == "chr6H"] - 199200000
delsnps$End[delsnps$chr == "chr7H"] <- delsnps$End[delsnps$chr == "chr7H"] - 254400000

## Plotting
delPlot <- ggplot(delsnps) +
	geom_segment(aes(x=delsnps$Start/1000000, y=1, xend=(delsnps$Start+5000000)/1000000, yend=1),size=16,color=delsnps$Colour) +
        xlim(0, 243.2) +
	facet_wrap(~ chr,nrow=1) +
	ylab("5")+
	theme_bw(base_size = baseSize) +
	theme(strip.background=element_rect(fill="white"), panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank())+
	theme(axis.title.x=element_blank(), axis.ticks=element_blank(), axis.text=element_blank()) +
	theme(strip.text.x=element_blank())+
	theme(legend.position = "none")

##===============================================
## Plotting Tajima's D
##===============================================
write("--------Plotting Tajima's D", stdout())
taj$chr <- factor(taj$CHROM,levels=chromosomeOrder)
taj$BIN_START[taj$chr == "chr1H"] <- taj$BIN_START[taj$chr == "chr1H"] - 151200000
taj$BIN_START[taj$chr == "chr2H"] <- taj$BIN_START[taj$chr == "chr2H"] - 256000000
taj$BIN_START[taj$chr == "chr3H"] <- taj$BIN_START[taj$chr == "chr3H"] - 224000000
taj$BIN_START[taj$chr == "chr4H"] <- taj$BIN_START[taj$chr == "chr4H"] - 163200000
taj$BIN_START[taj$chr == "chr5H"] <- taj$BIN_START[taj$chr == "chr5H"] - 168800000
taj$BIN_START[taj$chr == "chr6H"] <- taj$BIN_START[taj$chr == "chr6H"] - 199200000
taj$BIN_START[taj$chr == "chr7H"] <- taj$BIN_START[taj$chr == "chr7H"] - 254400000
dataTajSubset <- taj[grep ("chrUn", taj$chr, invert=TRUE), ]

rollAvgWindowSize <- as.integer(1000)

tajPlot<- ggplot(dataTajSubset, aes(BIN_START/1000, PI)) +
#geom_point(size = 0.3, colour="#d18884", aes(x=dataTajSubset$BIN_START/1000, y=dataTajSubset$spontaneums), na.rm=TRUE, alpha=0.8) +
#geom_point(size = 0.3, colour="#f0b299", aes(x=dataTajSubset$BIN_START/1000, y=dataTajSubset$landraces),   na.rm=TRUE, alpha=0.8) +
#geom_point(size = 0.3, colour="#7996b3", aes(x=dataTajSubset$BIN_START/1000, y=dataTajSubset$cultivars),   na.rm=TRUE, alpha=0.8) +
geom_line(size = lineWidth+1, colour="#b2182b", aes(x=dataTajSubset$BIN_START/1000, y=rollapply(as.numeric(as.character(dataTajSubset$spontaneums)), rollAvgWindowSize, mean, na.rm=TRUE, partial=TRUE, fill = NA)) ) +
geom_line(size = lineWidth+1, colour="#ef8a62", aes(x=dataTajSubset$BIN_START/1000, y=rollapply(as.numeric(as.character(dataTajSubset$landraces)), rollAvgWindowSize, mean, na.rm=TRUE, partial=TRUE, fill = NA)) ) +
geom_line(size = lineWidth+1, colour="#2166ac", aes(x=dataTajSubset$BIN_START/1000, y=rollapply(as.numeric(as.character(dataTajSubset$cultivars)), rollAvgWindowSize, mean, na.rm=TRUE, partial=TRUE, fill = NA)) ) +
geom_hline(aes(yintercept=1.19227))+
geom_hline(aes(yintercept=-2.04366))+
facet_wrap(~ dataTajSubset$chr,nrow=1) +
ylab("b")+
theme_bw(base_size = baseSize) +
theme(strip.background=element_blank(), panel.border=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank()) +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
theme(axis.title.y=element_text(angle = 0,face="bold"))+
theme(strip.text.x=element_blank())

##===============================================
## Plotting candidate genes
##===============================================
write("--------Plotting candidate genes",stdout())
genes$chr <- factor(genes$Chrom,levels=chromosomeOrder)

## Convert zone	3 coordinate
genes$Start[genes$chr == "chr1H"] <- genes$Start[genes$chr == "chr1H"] - 151200000
genes$Start[genes$chr == "chr2H"] <- genes$Start[genes$chr == "chr2H"] - 256000000
genes$Start[genes$chr == "chr3H"] <- genes$Start[genes$chr == "chr3H"] - 224000000
genes$Start[genes$chr == "chr4H"] <- genes$Start[genes$chr == "chr4H"] - 163200000
genes$Start[genes$chr == "chr5H"] <- genes$Start[genes$chr == "chr5H"] - 168800000
genes$Start[genes$chr == "chr6H"] <- genes$Start[genes$chr == "chr6H"] - 199200000
genes$Start[genes$chr == "chr7H"] <- genes$Start[genes$chr == "chr7H"] - 254400000
genes$End[genes$chr == "chr1H"] <- genes$End[genes$chr == "chr1H"] - 151200000
genes$End[genes$chr == "chr2H"] <- genes$End[genes$chr == "chr2H"] - 256000000
genes$End[genes$chr == "chr3H"] <- genes$End[genes$chr == "chr3H"] - 224000000
genes$End[genes$chr == "chr4H"] <- genes$End[genes$chr == "chr4H"] - 163200000
genes$End[genes$chr == "chr5H"] <- genes$End[genes$chr == "chr5H"] - 168800000
genes$End[genes$chr == "chr6H"] <- genes$End[genes$chr == "chr6H"] - 199200000
genes$End[genes$chr == "chr7H"] <- genes$End[genes$chr == "chr7H"] - 254400000

## Plotting
genePlot <- ggplot(genes) +
        geom_segment(aes(x=genes$Start/1000000, y=1, xend=(genes$Start+2500000)/1000000, yend=1),size=16,color=genes$Colour) +
        xlim(0, 243.2) +
        facet_wrap(~ chr,nrow=1) +
        ylab("b")+
        theme_bw(base_size = baseSize) +
        theme(strip.background=element_rect(fill="white"), panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank())+
        theme(axis.title.x=element_blank(), axis.ticks=element_blank(), axis.text=element_blank()) +
        theme(strip.text.x=element_blank())+
	theme(axis.title.y=element_text(angle = 0,face="bold"))+
        theme(legend.position = "none")

##===============================================
## Plotting the panel
##===============================================
write("--------Creating panel",stdout())
# Original figure contains all plot
#plotPanel <- plot_grid(coorPlot,haploPlot,muPlot1,muPlot2,muPlot3,delPlot,genePlot,align="v",nrow=7,rel_heights=c(0.3,1,1,1,1,0.6,0.6))
# New figure with just coordinate, haplotype count, and candidate genes
plotPanel <- plot_grid(coorPlot,delPlot,align="v",nrow=2,rel_heights=c(0.6,1))

write("--------Writing file to png",stdout())
png(outName,3000,500)
grid.draw(plotPanel)
dev.off()
write("--------Done",stdout())
