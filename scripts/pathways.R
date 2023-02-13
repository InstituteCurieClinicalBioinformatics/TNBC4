library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)
library(optparse)

pathwayEnrichment <- function(entrezID, cluster, outDir){
    pathways <- enrichGO(entrezID, 'org.Hs.eg.db', pvalueCutoff=1, qvalueCutoff=1)
    pathways = setReadable(pathways, 'org.Hs.eg.db')
    d = dotplot(pathways) + ggtitle(paste("Cluster ", cluster, sep=""))
    ggarrange(d, ncol=1) %>% ggexport(filename = paste(outDir, "/Cluster", cluster, ".jpg", sep=""), width = 2000, height = 1500, res = 200)
}

getEntrezID <- function(symbols){
    hs <- org.Hs.eg.db
    return(select(hs, keys = symbols, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL"))
}

getProcess <- function(processes, splitProcesses){
    splitProcesses = c(splitProcesses, as.list(strsplit(processes, '///')))
    return(splitProcesses)
}

getPathwayForOne <- function(data, group, outDir){
    entrezIDs = getEntrezID(rownames(data))$ENTREZID
    pathwayEnrichment(entrezIDs, group, outDir)
}

main <- function(inFile, outDir, threshold){
    dir.create(file.path(outDir, "pathwaysEnrichment"), showWarnings = FALSE)
    outDir = file.path(outDir, "pathwaysEnrichment")
    scaledData = read.csv(inFile, row.names=1)
    scaledData = scaledData[2:dim(scaledData)[2]]
    columnsPval = c("X1_VS_ALL", "X2_VS_ALL", "X3_VS_ALL", "X4_VS_ALL")
    for (i in 1:length(columnsPval)){
        if (dim(scaledData[scaledData[columnsPval[i]] < threshold, ])[1] == 0){
            message(paste("No significant differential expression observed for group ", i, " (pval threshold = ", threshold, ")", sep=""))
        }else{
            getPathwayForOne(scaledData[scaledData[columnsPval[i]] < threshold, ], i, outDir)
        }
    }
}

option_list = list(

make_option(c("-i", "--inFile"), type="character", help="Final file with pvalue"),
make_option(c("-o", "--outDir"), type="character", help="Output folder"),
make_option(c("-t", "--threshold"), type="double", help="Pvalue cutoff", default=0.05));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

main(opt$inFile, opt$outDir, opt$threshold)
