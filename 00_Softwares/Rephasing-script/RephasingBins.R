#!/usr/bin/env Rscript

options(stringsAsFactors=F)

args=commandArgs(T)

input <- args[1]
outdir <- args[2]

dat <- read.table(input,header=T)
dat <- as.matrix(dat)
missing_vec <- NULL
miss_threshold <- round(ncol(dat) * 0.3)
for (i in 1:nrow(dat)) {
	dat[i,] <- sub("A",0,dat[i,])
	dat[i,] <- sub("B",1,dat[i,])
	dat[i,] <- sub("U",NA,dat[i,])
	## missing rate, can be modified
	if (sum(is.na(dat[i,])) > miss_threshold ) {
		missing_vec <- c(missing_vec,i)
	}
}
if (length(missing_vec) != 0) {
	dat <- dat[-missing_vec,]
}

phasing_final_mat <- NULL
phasing_final_mat <- rbind(phasing_final_mat,dat[1,])
rownames(phasing_final_mat)[1] <- rownames(dat)[1]
dat <- dat[-1,]
list_count = 0
phasing_final_mat_list <- list()
wrong_phasing_bins <- NULL

while (nrow(dat) > 0) {
remove_vec <- NULL
mnames <- NULL
phasing_mat <- NULL
w_phasing_mat <- NULL
wmnames <- NULL

for (j in 1:nrow(phasing_final_mat)) {
	x <- phasing_final_mat[j,]
	for (i in 1:nrow(dat)) {
		vec <- x == dat[i,]
		ratio <- sum(vec,na.rm=T) / sum(!is.na(vec))
		if (ratio >= 0.92) {
			phasing_mat <- rbind(phasing_mat,dat[i,])
			mnames <- c(mnames,rownames(dat)[i])
			remove_vec <- c(remove_vec,i)
		}else if (ratio <= 0.08) {
			w_phasing_mat <- rbind(w_phasing_mat,dat[i,])
			wmnames <- c(wmnames,rownames(dat)[i])
			remove_vec <- c(remove_vec,i)
		}
	}
}

if (length(wmnames) == 0 & length(mnames) == 0) {
	list_count = list_count + 1
	phasing_final_mat_list[[list_count]] <- phasing_final_mat
	
	phasing_final_mat <- NULL
	phasing_final_mat <- rbind(phasing_final_mat,dat[1,])
	rownames(phasing_final_mat)[1] <- rownames(dat)[1]
	dat <- dat[-1,]
	
}else{
	
	if (length(wmnames) != 0) {
		w_phasing_mat <- ifelse(w_phasing_mat == 0,1,0)
		rownames(w_phasing_mat) <- wmnames
		rownames(phasing_mat) <- mnames
		w_phasing_mat <- w_phasing_mat[!duplicated(rownames(w_phasing_mat)),]
		phasing_mat <- phasing_mat[!duplicated(rownames(phasing_mat)),]
		
		wmnames <- wmnames[!duplicated(wmnames)]
		mnames <- mnames[!duplicated(mnames)]
		
		w_phasing_mat <- rbind(NULL,w_phasing_mat)
		phasing_mat <- rbind(NULL,phasing_mat)
		rownames(w_phasing_mat) <- wmnames
		rownames(phasing_mat) <- mnames
		
		phasing_final_mat <- rbind(phasing_final_mat,phasing_mat,w_phasing_mat)
		
		wrong_phasing_bins <- c(wrong_phasing_bins,wmnames)
	}else{
		rownames(phasing_mat) <- mnames
		phasing_mat <- phasing_mat[!duplicated(rownames(phasing_mat)),]
		
		mnames <- mnames[!duplicated(mnames)]
		phasing_mat <- rbind(NULL,phasing_mat)
		rownames(phasing_mat) <- mnames
		
		phasing_final_mat <- rbind(phasing_final_mat, phasing_mat)
	}
	remove_vec <- unique(remove_vec)
	if (nrow(dat) - length(remove_vec) == 1) {
		final_vec <- dat[-remove_vec,]
		final_m_name <- rownames(dat)[-remove_vec]
		final_gmt_name <- names(final_vec)
		dat <- matrix(final_vec,nrow=1)
		rownames(dat) <- final_m_name
		colnames(dat) <- final_gmt_name
	}else{
		dat <- dat[-remove_vec,]
	}
}
dat <- rbind(NULL,dat)
}

list_count = list_count + 1
phasing_final_mat_list[[list_count]] <- phasing_final_mat

final_res <- NULL

for (i in 1:length(phasing_final_mat_list)) {
	tmp <- phasing_final_mat_list[[i]]
	tmp <- ifelse(tmp == 0,"A","B")
	final_res <- rbind(final_res,tmp)
}

for (i in 1:nrow(final_res)) {
	final_res[i,][is.na(final_res[i,])] <- "U"
}


dim(final_res)

write.table(final_res, paste0(outdir, "/merged_binMarkers.re-phasing.co8.missing30.txt"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(wrong_phasing_bins, paste0(outdir, "/wrong_phasing_bins.txt"),row.names=F,col.names=F,quote=F,sep="\t")
