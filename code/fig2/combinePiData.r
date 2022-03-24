#Script for combining the pi data from population-wise output into a single file that can be used for plotting. 

library(tidyverse)
library(data.table)
library(zoo)

args<-commandArgs(TRUE)

message("read the input files")

cultivarData <- read.delim("pi_10k_windows.cultivars.windowed.pi", header=T)
landraceData <- read.delim("pi_10k_windows.landraces.windowed.pi", header=T)
wildBarleyData <- read.delim("pi_10k_windows.wildBarleys.windowed.pi", header=T)

#to merge the data, we need to create an additional column in each df that combines chromosome and bin start
cultivarData$id <- paste(cultivarData$CHROM, cultivarData$BIN_START, sep="_")
landraceData$id <- paste(landraceData$CHROM, landraceData$BIN_START, sep="_")
wildBarleyData$id <- paste(wildBarleyData$CHROM, wildBarleyData$BIN_START, sep="_")

str(cultivarData)
str(landraceData)
str(wildBarleyData)

colnames(cultivarData)

message("merge the data")
#merge the dataframes
#put all data frames into list
df_list <- list(cultivarData, landraceData, wildBarleyData)
#merge all data frames in list
data <- df_list %>% reduce(full_join, by='id')

#remove all rows where we have missing data
completeRowsOnly <- data[complete.cases(data), ]

message("sort the data")
dataSorted <- completeRowsOnly[with(completeRowsOnly, order(CHROM, BIN_START)), ]

#rename columns appropriately to make it backward compatible wirh existing plotting code
#CHROM   POS     cultivars       landraces       spontaneums
setnames(dataSorted, "BIN_START", "POS")
setnames(dataSorted, "PI.x", "cultivars")
setnames(dataSorted, "PI.y", "landraces")
setnames(dataSorted, "PI", "spontaneums")

#now pick only the columns we need
subset5cols <- data.frame(dataSorted$CHROM, dataSorted$POS, dataSorted$cultivars,dataSorted$landraces,dataSorted$spontaneums)
colnames(subset5cols) <- c("CHROM", "POS", "cultivars","landraces","spontaneums")
message("head(subset5cols)")
head(subset5cols)

#exclude chrUn
write("subset the data to exclude chrUn", stdout())
dataSubset <- subset5cols[grep ("chrUn", subset5cols$CHROM, invert=TRUE), ]

message("write data to tab delimited output file")
write.table(dataSubset,file="pi_10k_windows_combined.txt",row.names = F,col.names = T,quote = F,sep='\t')

message("workflow complete")

