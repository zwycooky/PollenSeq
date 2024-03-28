#!/usr/bin/env Rscript

library(Hapi)

args=commandArgs(T)

input1 <- args[1]
input2 <- args[2]


co_file <- read.table(input1,header=T)
gmt <- read.table(input2,header=T)

gmt <- gmt[gmt$confidence != "L" & gmt$rate > 0.8,]

hetDa <- gmt[,1:4]
ref <- hetDa$ref
alt <- hetDa$alt

#gmtDa <- gmt[, c(colnames(gmt)[grepl("[QSC]",colnames(gmt))], "hap1", "hap2") ]
gmtDa <- gmt[,5:(ncol(gmtDa)-3)]
gmtDa <- base2num(gmt = gmtDa, ref = ref, alt = alt)


merge_region <- function(dat) {
	
	if (ncol(dat) != 2) {
		stop("'dat' must be a 2 columns numeric matrix")
	}
	
	dat <- dat[order(dat[,1]),]
	
	merge_vec <- as.numeric(dat[-1,1] - dat[-nrow(dat),2] <= 0)
	if (sum(merge_vec) == 0) {
		return(dat)
	}else{
		start_site <- which(merge_vec[-1] - merge_vec[-length(merge_vec)] == 1) + 1
		end_site <- which(merge_vec[-length(merge_vec)] - merge_vec[-1] == 1) + 1
		
		if (merge_vec[1] == 1) {
			start_site <- c(1,start_site)
		}
		if (merge_vec[length(merge_vec)] == 1) {
			end_site <- c(end_site,length(merge_vec)+1)
		}
		
		merge_site <- cbind(start_site,end_site)
		
		dup_vec <- NULL
		merged_dat <- NULL
		for(i in 1:nrow(merge_site)) {
			tmp_vec <- merge_site[i,1]:merge_site[i,2]
			dup_vec <- c(dup_vec,tmp_vec)
			tmp_dat <- dat[tmp_vec,]
			new_s <- min(tmp_dat[,1])
			new_e <- max(tmp_dat[,2])
			merged_dat <- rbind(merged_dat,c(new_s,new_e))
		}
		
		dat <- dat[-dup_vec,]
		colnames(merged_dat) <- c("start","end")
		merged_dat <- rbind(dat,merged_dat)
		merged_dat <- merged_dat[order(merged_dat[,1]),]
		return(merged_dat)
	}
}

co_pos <- paste(co_file$start,co_file$end)
unique_co <- co_file[!duplicated(co_pos),]

merged_co <- merge_region(unique_co[,c("start","end")])
end_pos <- as.numeric(rownames(gmtDa)[nrow(gmtDa)])

binmap <- matrix(c(1,sort(unlist(merged_co)),end_pos),ncol=2,byrow=T)


## filter binmap ##
pos <- as.numeric(rownames(gmtDa))
remove_bin <- NULL
for (i in 1:nrow(binmap)) {
	bins <- binmap[i,1]
	bine <- binmap[i,2]
	bin_len <- bine - bins + 1
	bin_snp_num <- sum(pos > bins & pos <= bine)
	if ((bin_len < 200000 | bin_snp_num < 3) & i != 1 & i != nrow(binmap)) {
		remove_bin <- c(remove_bin, i)
	}
}

if (length(remove_bin) != 0) {
	binmap <- binmap[-remove_bin,]
}

if (length(nrow(binmap)) == 0 & length(binmap) != 0) {
	binmap <- rbind(binmap, NULL)
} else if (length(nrow(binmap)) == 0 & length(binmap) == 0) {
	stop(paste(input1,"too short for finding bins",sep=" "))
}




## get bin markers ##

phaseBin <- function (x) {
	x <- x[!is.na(x)]
	
	if (length(x) == 0) {
		return("U")
	}
	
	x <- abs(x)
	sum0 <- sum(x == 0)
	sum1 <- sum(x == 1)
	
	if (sum0 > sum1) {
		ratio <- sum0 / (sum1 + sum0)
		if (ratio > 0.85) {
			return("A")
		}else{
			return("U")
		}
	}else if (sum0 < sum1) {
		ratio <- sum1 / (sum1 + sum0)
		if (ratio > 0.85) {
			return("B")
		}else{
			return("U")
		}
	}else{
		return("U")
	}
	
}

pos <- as.numeric(rownames(gmtDa))
binMarkers <- matrix( 0, nrow=nrow(binmap), ncol=sum(grepl("[J,DASZ,QSC]",colnames(gmtDa))) )
for (i in 1:nrow(binmap)) {
	bins <- binmap[i,1]
	bine <- binmap[i,2]
	tmp_dat <- gmtDa[pos >= bins & pos <= bine,]
	tmp_gmt <- tmp_dat[,grepl("[J,DASZ,QSC]",colnames(tmp_dat))]
	hap <- tmp_dat$hap1
	
	phase_mat <- tmp_gmt - hap
	binMarkers[i,] <- apply(phase_mat,2,phaseBin)
}

mychr <- as.character(hetDa$chr[1])
Mname <- paste("m",mychr,1:nrow(binmap),sep="_")
rownames(binmap) <- Mname
colnames(binmap) <- c("start","end")
binmap <- data.frame(chr=rep(mychr,nrow(binmap)),binmap)

rownames(binMarkers) <- Mname
colnames(binMarkers) <- colnames(gmtDa)[grepl("[J,DASZ,QSC]",colnames(gmtDa))]

output_binmap <- paste(mychr,"binmap","txt",sep=".")
output_binMarkers <- paste(mychr,"binMarkers","txt",sep=".")

write.table(binmap,output_binmap,quote=F,row.names=T,col.names=T,sep="\t")
write.table(binMarkers,output_binMarkers,quote=F,row.names=T,col.names=T,sep="\t")
