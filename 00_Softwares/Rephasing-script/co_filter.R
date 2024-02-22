#!/usr/bin/env Rscript

options(stringsAsFactors=F)

Args <- commandArgs(T)
co_file <- Args[1]
reanchor_file <- Args[2]
output <- Args[3]

co <- read.table(co_file,header=F)
binmap <- read.table(reanchor_file,header=F)
rownames(binmap) <- binmap[,2] 

co <- cbind(co, binmap[co[,3],4], binmap[co[,4],3])
colnames(co) <- c("chromosome","gmt","bin1","bin2","start","end")

final_filter <- NULL
gmt <- unique(co$gmt)
for (i in 1:length(gmt)) {
	tmp <- co[co$gmt == gmt[i],]
	chr_vec <- unique(tmp$chromosome)
	for (j in 1:length(chr_vec)) {
		tmp1 <- tmp[tmp[,1] == chr_vec[j],]
		if (nrow(tmp1) == 1) {
			final_filter <- rbind(final_filter,tmp1)
		}else{
			vec <- as.numeric(tmp1$start[-1] - tmp1$end[-nrow(tmp1)] < 2000000)
			
			if (vec[1] == 1) {
				start_site <- c(1, which(vec[-1] - vec[-length(vec)] == 1) + 1 )
			}else{
				start_site <- which(vec[-1] - vec[-length(vec)] == 1) + 1
			}
			
			if (vec[length(vec)] == 1) { 
				end_site <- c(which(vec[-length(vec)] - vec[-1] == 1), length(vec)) + 1
			}else{
				end_site <- which(vec[-length(vec)] - vec[-1] == 1) + 1
			}
			
			remove_vec <- NULL
			mat <- cbind(start_site, end_site)
			if (nrow(mat) == 0) {
				final_filter <- rbind(final_filter,tmp1)
			}else{
			
				for (m in 1:nrow(mat)) {
					remove_vec <- c( remove_vec, -(mat[m,1] : mat[m,2]) )
				}
				tmp1 <- tmp1[remove_vec,]
				tmp2 <- tmp1
				final_filter <- rbind(final_filter,tmp1)
			}
		}
	}
}

write.table(final_filter, output,row.names=F,col.names=F,quote=F,sep="\t")
