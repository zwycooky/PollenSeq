#!/usr/bin/env Rscript

options(stringsAsFactors=F)
library(HMM)

Args <- commandArgs(T)
binmap_file <- Args[1]
Linkage <- Args[2]
output <- Args[3]

binmap <- read.table(binmap_file, header=T)
binmap[binmap == "U"] <- NA

gmtmiss <- 1:ncol(binmap)
for (i in 1:ncol(binmap)) {
        gmtmiss[i] <- sum(is.na(binmap[,i]))/nrow(binmap)
}
names(gmtmiss) <- colnames(binmap)

binmap <- binmap[,gmtmiss <= 0.1]

dat <- list()

file_names <- list.files(path = Linkage, pattern="lg")
for(i in 1:length(file_names)) {
        name <- gsub(".txt","",file_names[i])
        dat[[name]] <- read.table(paste0(Linkage,"/", file_names[i]),header=F,sep="\t")
}

co_detector <- function(gmt,hap1) {
        comp <- as.numeric(gmt == hap1)
        names(comp) <- 1:length(comp)
        comp <- comp[!is.na(comp)]


    hmm = initHMM(States=c("F","M"),
        Symbols=c("f","m"),
        transProbs=matrix(c(0.99999,0.00001,0.00001,0.99999),2),
        emissionProbs=matrix(c(0.99,0.01,0.01,0.99),2),
        startProbs = c(0.5,0.5))

        genoSymbol <- ifelse(comp==0,hmm$Symbols[1],hmm$Symbols[2])

    correctGeno <- viterbi(hmm, genoSymbol)
    correctGeno <- ifelse(correctGeno==hmm$States[1], 0, 1)
        names(correctGeno) <- names(comp)


        #lft <- which( abs(correctGeno[-1] - correctGeno[-length(correctGeno)]) == 1 )
        lft <- which( abs(comp[-1] - comp[-length(comp)]) == 1 )
        rht <- as.numeric(names(lft))
        lft <- as.numeric(names(comp[lft]))
        res <- cbind(lft,rht)
        return(res)
}

co_res <- NULL

for (i in 1:length(dat)) {
        tmp_dat <- dat[[i]]
        tmp_gmap <- binmap[tmp_dat[,1],]
        chr <- names(dat)[i]
        hap1 <- rep("A",nrow(dat[[i]]))
        for (j in 1:ncol(tmp_gmap)) {
                co_pos <- co_detector(tmp_gmap[,j],hap1)
                if (nrow(co_pos) != 0) {
                        gmt_id <- colnames(tmp_gmap)[j]
                        tmp_res <- cbind(rep(chr,nrow(co_pos)),rep(gmt_id,nrow(co_pos)),dat[[i]][,1][co_pos[,1]],dat[[i]][,1][co_pos[,2]])
                        co_res <- rbind(co_res,tmp_res)
                }
        }
}

## remove gmt which co numbers > 50 ##
print (names(which(table(co_res[,2]) > 50)))
co_res <- co_res[!co_res[,2] %in% names(which(table(co_res[,2]) > 50)),]
colnames(co_res) <- c("chr","sperm","left_bin_id","right_bin_id")

write.table(co_res, output,row.names=F,col.names=F,quote=F,sep="\t")
