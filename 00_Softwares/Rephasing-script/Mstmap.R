#!/usr/bin/env Rscript

library(ASMap)

options(stringsAsFactors=F)

args=commandArgs(T)

input <- args[1]
outdir <- args[2]

df <- read.csv(input, sep = "\t", header = TRUE)

dir.create(outdir)


results <- mstmap.data.frame(df,
                  pop.type = "DH",
                  dist.fun = "kosambi",
                  p.value = 0.00000001,
                  noMap.dist = 15.0,
                  noMap.size = 0,
                  miss.thresh = 1.00,
                  objective.fun = "COUNT",
                  detectBadData = TRUE,trace = TRUE
                  )


for (group_name in names(results$geno)) {
  map_data <- round(results$geno[[group_name]]$map, 3)
  file_name <- paste0("lg", gsub("L", "", group_name), ".txt")  # Constructing the file name
  full_file_path <- file.path(outdir, file_name)  # Constructing the full file path
  write.table(map_data, file = full_file_path,  sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
}



png(paste0(outdir,"/MStmap_result_sumamry.png"), width = 15, height = 10)
plot(results)
dev.off()

png(paste0(outdir, "/MStmap_Missing.png"), width = 6, height = 6)
plotMissing(results)
dev.off()

png(paste0(outdir, "/MStmap_GeneticMap.png"),  width = 6, height = 6)
plotMap(results)
dev.off()

png(paste0(outdir, "/MStmap_pheno.png"),  width = 6, height = 6)
plotPheno(results, pheno.col=1)
dev.off()
