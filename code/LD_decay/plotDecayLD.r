args<-commandArgs(TRUE)

write("## load packages", stdout())
library(zoo)
library(ggplot2)

write("## read the input file", stdout())
ld <- read.table(args[1], header=T,sep="\t")

write("## output to png file", stdout())

ldPlot <- ggplot(ld,aes(x=DIST/1000,y=R2,colour=GROUP))+
geom_smooth(method="auto",size=3) +
labs(x="Distance (Kbp)",y=expression(Average~ r^{2}))+
theme_gray(base_size=60)+
theme(panel.border=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank()) +
theme(legend.position = "none")+
scale_color_manual(values = c("#4575b4", "#d73027" , "#fc8d59"))+

png(args[2],2000,1000)
print(ldPlot)
dev.off()


write("done", stdout())
