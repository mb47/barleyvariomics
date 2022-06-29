args<-commandArgs(TRUE)

write("read the input file", stdout())
data <- read.table(args[1],sep="\t",header=TRUE)
outName <- args[2]

write("get packages", stdout())
library(ggplot2)
library(rworldmap)
library(ggmap)

write("prepare palette", stdout())
palette(c("#2166ac","#67a9cf","#ef8a62","#fddbc7","#b2182b"))
# pallette description: dark blue, blue, dark orange, orange, red
# for:			c6row, c2row, l6row, l2row, wild
## somehow the order in legend is different from that graph itself thus need manual adjustment (i.e. colour in legend 2 and 6 rows is inverted)

write("plotting", stdout())
# get base map
newmap <- getMap(resolution = "li")

# plotting
groupID <- unique(data$GROUP.ROW)
png(outName,6000,3000)
# Option 1: all lands in white
#mapCountryData(newmap,xlim=c(min(data$LONGTITUDE)-20,max(data$LONGTITUDE)+20),ylim=c(min(data$LATTITUDE)-20,max(data$LATTITUDE)+20),oceanCol='azure',colourPalette=c('white','white','white','white','white','white','white'),missingCountryCol='white',mapTitle='',addLegend=FALSE)
# Option 2: all lands in brownish (lemonchiffon2)
#mapCountryData(newmap,xlim=c(min(data$LONGTITUDE)-20,max(data$LONGTITUDE)+20),ylim=c(min(data$LATTITUDE)-20,max(data$LATTITUDE)+20),oceanCol='azure',colourPalette=c('lemonchiffon2','lemonchiffon2','lemonchiffon2','lemonchiffon2','lemonchiffon2','lemonchiffon2','lemonchiffon2'),missingCountryCol='lemonchiffon2',mapTitle='',addLegend=FALSE)
# Option 3: all lands in gray (gray), without extension of boundary
mapCountryData(newmap,xlim=c(min(data$LONGTITUDE),max(data$LONGTITUDE)),ylim=c(min(data$LATTITUDE),max(data$LATTITUDE)),oceanCol='white',colourPalette=c('gray','gray','gray','gray','gray','gray','gray'),missingCountryCol='gray',mapTitle='',addLegend=FALSE)
## with jitter
points(jitter(data$LONGTITUDE,40),jitter(data$LATTITUDE,40),cex=5,pch=1,col=data$COLOUR,lwd=8)
## w/o jitter
#points(data$LONGTITUDE,data$LATTITUDE,cex=5,pch=1,col=data$COLOUR,lwd=8)
# uncomment to add legend
legend("bottomleft",pch=19,col=groupID,legend=groupID,cex=6,pt.cex=6)
dev.off()
