{
    blockPoint <- which(draftHap$cvlink >= cvlink)
    if (length(blockPoint) == 0) {
        message("No region is required for proofreading !")
        hap <- draftHap$hap
        names(hap) <- rownames(draftHap)
        return(hap)
    }
    else {
        start <- 1
        hapBlock <- list()
        for (i in 1:length(blockPoint)) {
            end <- blockPoint[i] - 1
            hapBlock[[i]] <- draftHap[start:end, ]
            start <- blockPoint[i]
            if (i == length(blockPoint)) {
                end <- nrow(draftHap)
                hapBlock[[i + 1]] <- draftHap[start:end, ]
            }
        }
        message(paste("Size of blocks: ", paste(unlist(lapply(hapBlock, 
            nrow)), collapse = ","), "\n", sep = ""))
        filter <- which(lapply(hapBlock, nrow) < smallBlock)
        message(paste(length(filter), " blocks are removed !\n", 
            sep = ""))
        hapBlock[filter] <- NULL
        gmt <- gmtFrame[rownames(draftHap), ]
        currentHap <- hapBlock[[1]]
		
		if (length(hapBlock) == 1) {
			hap <- currentHap$hap
			names(hap) <- rownames(currentHap)
			
			filter_hap7 <- which(hap==7)
			if (length(filter_hap7) > 0) {
				hap <- hap[-filter_hap7,]
			}
			
			return(hap)
		}else{
			for (i in 2:(length(hapBlock))) {
				currentHap <- MPRFun2(gmt, currentHap, hapBlock[[i]])
			}
		}
		
        hap <- currentHap$hap
        names(hap) <- rownames(currentHap)
        return(hap)
    }
}

MPRFun2 <- function(gmt, hapBlock1, hapBlock2, nSNP=100) {
    
    filter <- which(hapBlock1$hap==7)
    if (length(filter) > 0) {
        hapBlock1 <- hapBlock1[-filter,]
    }
    
    
    filter <- which(hapBlock2$hap==7)
    if (length(filter) > 0) {
        hapBlock2 <- hapBlock2[-filter,]
    }
    
    
    if (nrow(hapBlock1)>nSNP) {
        sites <- rownames(hapBlock1)[(nrow(hapBlock1)-nSNP+1):nrow(hapBlock1)]
        raw1 <- gmt[sites,]
        hap1 <- hapBlock1[(nrow(hapBlock1)-nSNP+1):nrow(hapBlock1),]
        
    } else {
        raw1 <- gmt[rownames(hapBlock1),]
        hap1 <- hapBlock1
    }
    
    if (nrow(hapBlock2)>nSNP) {
        sites <- rownames(hapBlock2)[1:nSNP]
        raw2 <- gmt[sites,]
        hap2 <- hapBlock2[1:nSNP,]
    } else {
        raw2 <- gmt[rownames(hapBlock2),]
        hap2 <- hapBlock2
    }
    
    raw <- rbind(raw1, raw2)
    
    hap11 <- rbind(hap1, hap2)
    
    hap2$hap <- flipFun2(hap2$hap)
    hap22 <- rbind(hap1, hap2)
    
    sum11 <- sum(apply(raw, 2, function(v) cvCountFun2(hap11$hap, v)))
    sum22 <- sum(apply(raw, 2, function(v) cvCountFun2(hap22$hap, v)))
    
    if (sum11 > sum22) {
        hapBlock2$hap <- flipFun2(hapBlock2$hap)
    }
    
    hapBlock <- rbind(hapBlock1, hapBlock2)
    
    message (paste('Number of crossovers given haplotype 1/2: ', 
        sum11, '/', sum22, sep=''))
    
    return (hapBlock)
}

flipFun2 <- function(v){
    v2 <- ifelse(v==7,7,ifelse(v==0, 1, 0))
    return (v2)
}

cvCountFun2 <- function(gmt1, gmt2) {
    
    idComp <- gmt1 == gmt2
    idComp <- as.numeric(idComp[!is.na(idComp)])
    
    if (length(idComp)<=1) {
        return (0)
    }
    
    v1 <- idComp[-1]
    v2 <- idComp[-length(idComp)]
    vd <- v1-v2
    
    cv <- sum(vd != 0)
    return (cv)
}

