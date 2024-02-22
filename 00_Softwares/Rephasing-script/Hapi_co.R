#!/usr/bin/env Rscript

library('Hapi')
library('HMM')

args=commandArgs(T)

input <- args[1]
output <- args[2]

hapdat <- read.table(input,header=T, colClasses = c ("character"))
hapdat <- hapdat[hapdat$confidence != "L" & hapdat$rate > 0.8,]

hap <- hapdat[,c("hap1","hap2")]
gmt <- hapdat[,!colnames(hapdat) %in% c("chr", "pos", "ref", "alt", "hap1", "hap2", "total", "rate", "confidence")]


cvOutput <- hapiIdentifyCV(hap = hap, gmt = gmt)

write.table(cvOutput,output,sep="\t",col.names=T,row.names=F,quote=F)
