#Plot for making violin plots of on-/off-target stats

#input data looks like this:
# onTarget	25069
# onTarget	3470
# onTarget	17042
# offTarget	210
# offTarget	5682
# offTarget	12445

args<-commandArgs(TRUE)

library(ggplot2)
library(ggpubr)
library(scales)

#read the data
data <- read.table(args[1], header = FALSE)

#PNG output file
png(args[2], type="cairo")

#violin plot by on-/off-target category column
ggplot(data, aes(x=data$V1, y=data$V2, fill=data$V1)) + 
geom_violin() + 
scale_fill_manual(values=c("#999999", "#999999")) + 
geom_boxplot(width=0.1) +
scale_y_continuous(trans = "log10", breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)), name = args[3]) +
xlab("") + 
theme(text = element_text(size=20), legend.position="none")

dev.off()

