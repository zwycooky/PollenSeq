#!/usr/bin/env Rscript

library(HMM)
options(stringsAsFactors=F)

args=commandArgs(T)
input1 <- args[1]
input2 <- args[2]

dat <- read.table(input1,header=T)

## process chr length file ##
chr_lens <- read.table(input2,header=F)
rownames(chr_lens) <- chr_lens[,1]
## ----------------------- ##

mycontig <- as.character(dat[1,1])

haps <- dat[,c("hap1","hap2")]
#hapdat <- dat[,grepl("[QSC]",colnames(dat))]
hapdat <- dat[,5:(ncol(dat)-5)]

## co == 0 ##
phasing_dat <- matrix(0,nrow=1,ncol=ncol(hapdat))
for (i in 1:ncol(hapdat)) {
	comp <- as.numeric(hapdat[,i] == haps[,1])
	comp <- comp[!is.na(comp)]
	if (length(comp) == 0) {
		phasing_dat[1,i] <- 'U'
	}else if ( sum(comp == 1)/length(comp) > 0.85 ) {
		phasing_dat[1,i] <- 'A'
	}else if ( sum(comp == 0)/length(comp) > 0.85 ) {
		phasing_dat[1,i] <- 'B'
	}else{
		phasing_dat[1,i] <- 'U'
	}
}

colnames(phasing_dat) <- colnames(hapdat)
rownames(phasing_dat) <- paste("m",mycontig,1,sep="_")

binMarkers_output <- paste(mycontig,"binMarkers","txt",sep=".")
write.table(phasing_dat, binMarkers_output, col.names=T, row.names=T, quote=F, sep="\t")

output_binmap <- paste(mycontig,"binmap","txt",sep=".")
bin_len <- as.numeric(chr_lens[mycontig,2])
marker_id <- paste("m",mycontig,1, sep="_")
binmap <- cbind(mycontig, 1, bin_len)
colnames(binmap) <- c("chr","start","end")
rownames(binmap) <- marker_id
write.table(binmap,output_binmap,quote=F,row.names=T,col.names=T,sep="\t")
