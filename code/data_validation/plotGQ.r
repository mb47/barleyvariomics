args<-commandArgs(TRUE)

library(ggplot2)
library(scales)

outName <- args[2]

data  <- read.table(args[1],sep="",header=TRUE,row.names=NULL)
write("Read in data completed", stdout())

p <- ggplot(data, aes(x=GQ)) +
geom_histogram(fill="#999999",colour="black",lwd=1) +
facet_wrap(~data$TYPE) +
scale_y_continuous(trans = "log10", breaks = trans_breaks("log10", function(x) 10^x), limits=c(10^0,10^9),labels = trans_format("log10", math_format(10^.x))) +
#scale_y_continuous(labels = scales::comma) +
xlab("GQ")+
ylab("Count (log)")+
theme_gray(base_size=40)

write("Write to png file", stdout())
png(outName,2000,1000)
print(p)
dev.off()
