args<-commandArgs(TRUE)

## I/O
write("--------Read the input file", stdout())
data <- read.csv(args[1],header=T)
# Define chromosomes to iterate
chroms <- c("chr1H","chr2H","chr3H","chr4H","chr5H","chr6H","chr7H")
# Define groups to iterate
groups <- c("Cultivar","Landrace","Spontaneum")
# Define colour pallete
Pal <- as.character("Spectral")

## Libraries
write("--------Get packages", stdout())
library(ggplot2)
library(rworldmap)
library(ggmap)
library(RColorBrewer)
library(maptools)
library(scatterpie)
library(ggthemes)
library(mapproj)
library(cowplot)
library(gridGraphics)

## get base map
newmap <- getMap(resolution = "li")

# plotting
for (chrom in chroms)
{
	write(paste("--------Plotting",chrom), stdout())
	write("----------Plot all groups",stdout())
	
	#### Get color pallete
	## Get maximum group number
	maxGroupNo <- max(get(chrom,data),na.rm=T)
	## Get colour pallet
	#palette(colorRampPalette(brewer.pal(n=11,name=Pal))(maxGroupNo))
	palette(brewer.pal(n=maxGroupNo,name=Pal))
	## PieChart: Get pallet information for pie plot
	Pal_pie <- palette(colorRampPalette(brewer.pal(n=11,name=Pal))(maxGroupNo))
	## PieChart: Get list of groups
	listTemp <- 1:maxGroupNo
	list <- paste("Group",listTemp,sep="")

	#### Plot histogram
	outName <- paste(chrom,"_histogram.png",sep="")
	p <- ggplot(data,aes_string(x=chrom,fill=chrom))+
		geom_bar(stat="count")+
		labs(x="Groups",y="Number of individuals")+
		scale_x_continuous(limits=c(0,maxGroupNo+1),breaks=1:maxGroupNo,labels=1:maxGroupNo)+
		theme_gray(base_size=60)+
		theme(axis.text.x = element_text(angle=-45))+
		theme(panel.border=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank()) +
		theme(legend.position = "none")
	png(outName,2000,2000)
	print(p)
	dev.off()
	

	## Start plotting
	outName <- paste(chrom,"_allGroups.png",sep="")
        png(outName,6000,3000)

	## Option 3: all lands in gray (gray), without extension of boundary
	title=''
	mapCountryData(newmap,xlim=c(min(data$LONGTITUDE),max(data$LONGTITUDE)),ylim=c(min(data$LATTITUDE),max(data$LATTITUDE)),oceanCol='white',colourPalette=c('gray','gray','gray','gray','gray','gray','gray'),missingCountryCol='gray',mapTitle=title,addLegend=FALSE)
	## Add data point with jitter
	points(jitter(data$LONGTITUDE,40),jitter(data$LATTITUDE,40),cex=5,col=get(chrom,data),lwd=8,pch=get(chrom,data))
	
	#### Legend
	## Get group ID	(for legend)
        groupID <- sort(unique(get(chrom,data)))
	## Add legend
	legend("bottomleft",pch=groupID,col=groupID,legend=groupID,cex=6,pt.cex=6,lwd=4)
	
	dev.off()

	#### Plotting different origins
	for (group in groups)
	{
		write(paste("----------Plot",group),stdout())
		worldmap <- map_data("world")

		######################
		#### Scatter plot ####
		######################
		#### Subset data
		#dataSubset <- data[grep(group,data$CATEGORY),]
                ## Calculate xlim and assign fix radius
                #minx <- min(dataSubset$LONGTITUDE)
                #maxx <- max(dataSubset$LONGTITUDE)
		#miny <- min(dataSubset$LATTITUDE)
		#maxy <- max(dataSubset$LATTITUDE)

                #mapplot1 <- ggplot(worldmap,aes(long,lat))+
                #        geom_map(map=worldmap,aes(map_id=region),col="black",fill="gray",size=1)+
                #        geom_point(data=dataSubset,aes(x=jitter(LONGTITUDE,40),y=jitter(LATTITUDE,40),color=get(chrom,dataSubset)),shape=get(chrom,dataSubset))+
                #        #scale_color_manual(values=Pal_pie)+
                #        coord_equal(xlim=c(minx-10,maxx+10),ylim=c(miny-10,maxy+10))+
                #        theme_bw(base_size=60)+
                #        theme(panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank())+
                #        theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text=element_blank())
		#tempOut=paste(chrom,group,"_test.png",sep="")
		#png(tempOut,1000,500)
		#print(mapplot1)
		#dev.off()
		
       		## Get maximum group number
        	#maxGroupNo <- max(get(chrom,dataSubset),na.rm=T)
        	## Get colour pallet
        	#palette(colorRampPalette(brewer.pal(n=11,name=Pal))(maxGroupNo))
		#palette <- palette(colorRampPalette(brewer.pal(n=11,name=Pal))(maxGroupNo))
 	 	## Start plotting
        	#outName <- paste(chrom,"_",group,".png",sep="")
        	#png(outName,1500,500)
		
		## Option 3: all lands in gray (gray), without extension of boundary
		#title=paste(chrom,group,sep="-")
	        #mapCountryData(newmap,xlim=c(min(dataSubset$LONGTITUDE),max(dataSubset$LONGTITUDE)),ylim=c(min(dataSubset$LATTITUDE),max(dataSubset$LATTITUDE)),oceanCol='white',colourPalette=c('gray','gray','gray','gray','gray','gray','gray','gray'),missingCountryCol='gray',mapTitle=title,addLegend=FALSE)
        	## Add data point with jitter
        	#points(jitter(dataSubset$LONGTITUDE,40),jitter(dataSubset$LATTITUDE,40),cex=2,col=get(chrom,dataSubset),lwd=2,pch=get(chrom,dataSubset))
		
	        #### Legend
        	## Get group ID (for legend)
        	#groupID <- sort(unique(get(chrom,dataSubset)))
        	## Add legend
        	#legend("bottomleft",pch=groupID,col=groupID,legend=groupID,cex=1,pt.cex=1,lwd=1)
		
		##################
		#### Pie plot ####
		##################
		fileToRead <- paste("dataframe/df_",chrom,".",group,".tsv",sep="")
		pieData <- read.table(fileToRead,sep="\t",header=T,quote="",fill=FALSE)
		worldmap <- map_data("world")
		## Calculate xlim and assign fix radius
		if (group != "Spontaneum"){
			minx <- min(data$LONGTITUDE)
			maxx <- max(data$LONGTITUDE)
			miny <- min(data$LATTITUDE)
			maxy <- max(data$LATTITUDE)
		} else {
                        minx <- min(pieData$Long)
                        maxx <- max(pieData$Long)
                        miny <- min(pieData$Lat)
                        maxy <- max(pieData$Lat)
		}
		rad <- as.integer(3)
		## Plotting		
		## For representing population size as radius, use scatterpie(aes(r=sqrt(Radius/3.14)))
		#mapplot2 <- ggplot(worldmap,aes(long,lat))+
                #        geom_map(map=worldmap,aes(map_id=region),col="black",fill="gray",size=1)+
                #        geom_scatterpie(aes(x=Long,y=Lat,group=Country,r=rad),data=pieData,cols=colnames(pieData[,5:ncol(pieData)]),lwd=1)+
                #        scale_fill_manual(breaks=list,values=Pal_pie)+
		#	coord_equal(xlim=c(minx-10,maxx+10),ylim=c(miny-10,maxy+10))+
                #        theme_bw(base_size=60)+
                #        theme(panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank())+
                #        theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text=element_blank())
		
		#outName <- paste(chrom,"_",group,".png",sep="")
		#png(outName,3000,1000)
		#print(mapplot2)
		#dev.off()

		#### Plot panel
		#plotPanel <- plot_grid(mapplot1,mapplot2,nrow=1,ncol=2)
		#outName <- paste(chrom,"_",group,".png",sep="")
		#png(outName,3000,1000)
		#grid.draw(plotPanel)
		#dev.off()
	}
	write("------------Done",stdout())
}

write("------------------------",stdout())
write("--------Warning messages",stdout())
warnings()
