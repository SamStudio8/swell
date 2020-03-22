library(tidyverse)
library(ggplot2)

df <- read.table("table", head=F, sep=' ', na=c("NaN"))

df$median_depth <- df$V7
df$is_20x <- df$V7 >= 20
df[df$median_depth == 0,]$median_depth <- 1
bed <- read.table("nCoV-2019.scheme.bed.left", head=F, sep='\t')

p <- ggplot(df, aes(x=V2, y=V1, fill=log10(median_depth))) + geom_tile() + theme(text=element_text(size=20)) + theme(legend.position="right") + scale_x_continuous(expand=c(0,0), breaks=seq(1,98), labels=c(paste(seq(1,98), '-', bed$V3))) + theme(axis.text.x=element_text(angle=90, hjust=1)) + scale_fill_gradientn(colours = rainbow(5))

q <- ggplot(df, aes(x=V2, y=V1, colour="black", fill=is_20x)) + geom_tile() + theme(text=element_text(size=20)) + theme(legend.position="right") + scale_x_continuous(expand=c(0,0), breaks=seq(1,98), labels=c(paste(seq(1,98), '-', bed$V3))) + theme(axis.text.x=element_text(angle=90, hjust=1))
