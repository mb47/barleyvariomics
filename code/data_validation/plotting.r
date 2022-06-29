data  <- read.table("othersMarkers.GQ.GQ.FORMAT",sep="\t",header=TRUE,na.strings=".")
data2 <- read.table("lowMAF.GQ.GQ.FORMAT",sep="\t",header=TRUE,na.strings=".")

data$AVG <- rowMeans(data[,3:817])
data2$AVG <- rowMeans(data2[,3:817])

data$group <- paste("Other SNPs")
data2$group <- paste("Low MAF SNPs")

library(ggplot2)

png("GQ_boxplot.png",2000,2000)
p <- ggplot()+
geom_violin(fill="#999999",aes(x=data$group,y=data$AVG),lwd=2)+
geom_violin(fill="#999999",aes(x=data2$group,y=data2$AVG),lwd=2)+
geom_boxplot(aes(x=data$group,y=data$AVG),width=0.1,fill="#999999",size=2)+
geom_boxplot(aes(x=data2$group,y=data2$AVG),width=0.1,fill="#999999",size=2)+
xlab("")+
ylab("Average GQ score")+
theme_gray(base_size=80)

print(p)
dev.off()
