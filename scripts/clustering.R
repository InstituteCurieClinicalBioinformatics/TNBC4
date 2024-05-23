library(ggplot2)
library(ggpubr)
library(FactoMineR)
library(factoextra)
library(dplyr)
library(gridExtra)
require(grid)
library(biomaRt)
library(affy)
library(limma)
library(stringr)
library(rstatix)
library(optparse)

formateForGGplot <- function(samples, graph){
    sampleNames = c()
    values = c()
    for (i in 1:length(samples)){
        sampleNames = c(sampleNames,rep(samples[i],dim(graph)[1]))
        values = c(values, graph[,i])
    }
    df = data.frame(Intensity=values,Sample=sampleNames)
    return(df)
}


#Plusieurs probes pour un meme gene. Recherche des lignes pour un meme gene et calcul de la moyenne pour ce gene
mergeIDvalue <- function(symbol, df, samples, genesProcessed){
    if (!(is.na(symbol)) & (!(symbol %in% genesProcessed))){
        sub = subset(df, Gene.symbol == symbol)
        duplicateMeans = as.data.frame(t(as.data.frame(colMeans(subset(sub, select = samples)))))
        rownames(duplicateMeans) = symbol
        return(duplicateMeans)
    }
}

#lapply genere un out of memory si lance sur toute la liste de genes. Donc on split la liste de genes puis on lance lapply sur les parties splitees
splitDF <- function(df, factor, samples){
    symbolWithoutDup = unique(df$Gene.symbol)
    genesProcessed = c()
    dfWithoutDup = list()
    len = round(length(symbolWithoutDup) / factor)
    for (i in 1:factor){
        start = len*(i-1)
        if (((length(symbolWithoutDup) %% factor) != 0) & (i == factor)){
            end = length(symbolWithoutDup)
        }else{
            end = len*i
        }
        symbols = symbolWithoutDup[start:end]
        duplicateMeans = lapply(symbols, mergeIDvalue, df = df, samples = samples, genesProcessed = genesProcessed)
        dfWithoutDup[[i]] = do.call("rbind", duplicateMeans)
        genesProcessed = c(genesProcessed, symbols)
    }
    return(dfWithoutDup)
}

findBestKmeans <- function(df, k, iter){
    lowestWSS = 100000000000000000000000000000000000000000
    keepKM = NULL
    for (i in 1:1000){
        res.km = kmeans(scale(as.data.frame(t(df))), k, iter)
        if (res.km$tot.withinss < lowestWSS){
            lowestWSS = res.km$tot.withinss
            keepKM = res.km
        }
    }
    return(keepKM)
}

mergeIndividual <- function(numHisto, df, data){
    test = as.data.frame(subset(df, Num_Histo == numHisto))
    if ((test$Batch2 != "") & (test$Batch2 != 0)){
        sampleID1 = unlist(strsplit(test$Batch1, "_"))
        sampleID2 = unlist(strsplit(test$Batch2, "_"))
        toMerge = c(paste(sampleID1[1], sampleID1[2], sep="_"), paste(sampleID2[1], sampleID2[2], sep="_"))
        individualMerged = as.data.frame(t(as.data.frame(colMeans(subset(data, rownames(data) %in% toMerge)))))
        rownames(individualMerged) = toMerge[1]
        return(individualMerged)
    }
}

mergeDf <- function(dfList){
    count = 1
    df = dfList[[count]]
    while (count < length(dfList)){
        count = count + 1
        merged = merge(df, dfList[[count]], by=0)
        if (dim(merged)[1] != 0){
            rownames(merged) = merged$Row.names
            df = dplyr::select(merged, -c('Row.names'))
        }
    }
    return(df)
}

toRemove <- function(sampleIDs){
    remove = c()
    for (sample in sampleIDs){
        id = unlist(strsplit(sample, "_"))
        remove = c(remove, paste(id[1], id[2], sep="_"))
    }
    return(remove)
}

meanBoxplot <- function(data, genes, group, samplesID, km, conclusion){
    subGroup = subset(data, rownames(data) %in% genes)
    if (dim(subGroup)[1] == 0){
        message(paste("Can't boxplot info for ", group, " genes\nGenes are not available in data.", sep=""))
    }else{
        geneVal = as.data.frame(colMeans(subGroup))
        colnames(geneVal) = c('mean')
        geneVal = as.data.frame(t(geneVal))
        geneVal = formateForGGplot(samplesID, geneVal)
        geneVal$cluster = km$cluster
        my_comparisons = list(c(1,2), c(1,3), c(1,4), c(2,3), c(2,4), c(3,4))
        p = compare_means(Intensity ~ cluster, geneVal, method = "kruskal.test", paired=FALSE)
        if (group %in% c("IM", "LAR", "MES", "BLIS")){
            conclusion = c(conclusion, attributeSubtype(geneVal, p$p.adj<0.05, group))
        }
        if (p$p.adj < 0.05){
            genePlot <- ggboxplot(geneVal, x = "cluster", y = "Intensity", add = "jitter", palette = "jco", color = "cluster", main = group)+ stat_compare_means(label.y = round(max(geneVal$Intensity)) + 5) + stat_compare_means(label = "p.signif", comparisons = my_comparisons, hide.ns = TRUE) + rremove("legend")
        }
        else{
            genePlot <- ggboxplot(geneVal, x = "cluster", y = "Intensity", add = "jitter", palette = "jco", color = "cluster", main = group)+ stat_compare_means(label.y = round(max(geneVal$Intensity)) + 5) + rremove("legend")
        }
        return(list("genePlot" = genePlot, "conclusion" = conclusion))
    }
}

attributeSubtype <- function(geneVal, significative, subtype){
    geneVal$cluster = as.factor(geneVal$cluster)
    medianCluster1 = median(subset(geneVal, cluster == 1)$Intensity)
    medianCluster2 = median(subset(geneVal, cluster == 2)$Intensity)
    medianCluster3 = median(subset(geneVal, cluster == 3)$Intensity)
    medianCluster4 = median(subset(geneVal, cluster == 4)$Intensity)
    higherCluster = which.max(c(medianCluster1, medianCluster2, medianCluster3, medianCluster4))
    comparaisons = as.data.frame(tukey_hsd(geneVal, Intensity ~ cluster))
    comparaisons$groupes = paste(comparaisons$group1, comparaisons$group2, sep = "_")
    comparaisons = subset(comparaisons, grepl(higherCluster, groupes))
    if (significative){
        comparaisons = as.data.frame(tukey_hsd(geneVal, Intensity ~ cluster))
        comparaisons$groupes = paste(comparaisons$group1, comparaisons$group2, sep = "_")
        comparaisons = subset(comparaisons, grepl(higherCluster, groupes))
        if (length(comparaisons$p.adj[comparaisons$p.adj < 0.05]) == 3){
            return(paste(higherCluster, "=", subtype, sep = ""))
        }
        else if (length(comparaisons$p.adj[comparaisons$p.adj < 0.05]) == 0){
            return(paste(higherCluster, "=", subtype, "(ns)", sep = ""))
        }
        else{
            return(paste(higherCluster, "=", subtype, "(nas)", sep = ""))
        }
    }
    else{
        return(paste(higherCluster, "=", subtype, "(ns)", sep = ""))
    }
}

getData <- function(files, samplesID){
    data = oligo::read.celfiles(files)
    colnames(data) = samplesID
    return(data)
}

computeQC <- function(data, outDir, cohorteName, files){
    fitData = fitProbeLevelModel(data)
    dir.create(file.path(outDir, "QC"), showWarnings = FALSE)
    QCdir = file.path(outDir, "QC")

    ###RLE
    filenameRLE = file.path(QCdir, paste(cohorteName, "RLE.jpg", sep="_"))
    rle = oligo::RLE(fitData, "values")
    df = formateForGGplot(colnames(data), rle)
    RLEbox = ggboxplot(df, x = "Sample", y = "Intensity", add = "jitter", palette = get_palette("jco", length(colnames(data))), color = "Sample", size=0.5) + rremove("legend") + rotate_x_text() +  font("xy.text", size = 10)
    ggarrange(RLEbox, ncol=1) %>% ggexport(filename = filenameRLE, width=668, height=662, res = 120)

    ###NUSE
    filenameNUSE = file.path(QCdir, paste(cohorteName, "NUSE.jpg", sep="_"))
    nuse = oligo::NUSE(fitData, "values")
    df = formateForGGplot(colnames(data), nuse)
    NUSEbox = ggboxplot(df, x = "Sample", y = "Intensity", add = "jitter", palette = get_palette("jco", length(colnames(data))), color = "Sample", size=0.5) + rremove("legend") + rotate_x_text() +  font("xy.text", size = 10)
    ggarrange(NUSEbox, ncol=1) %>% ggexport(filename = filenameNUSE, width=668, height=662, res = 120)
}

preProcessing <- function(data, chipAnnotation, nbSplit){
    samplesID = rownames(data@phenoData)
    data.rma = oligo::rma(data)
    data.matrix = exprs(data.rma)#[1:10,]
    ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
    message("Retrieving annotation")
    data.annotation = getBM(attributes = c(chipAnnotation, "hgnc_symbol"), mart = ensembl)
    data.annotation = data.annotation[data.annotation[, 1] != "", ] #remove line with no chip annotation
    data.annotation = data.annotation[data.annotation[, 2] != "", ] #remove line with no hgnc symbol
    colnames(data.annotation) = c("ID", "Gene.symbol")

    ###Add annotation
    df = as.data.frame(data.matrix)
    df$ID = rownames(df)
    df$ID = str_split_fixed(df$ID, "\\.1", 2)[,1]
    data.annotated = merge(df, data.annotation, by=c("ID"))
    finalData = dplyr::select(data.annotated, samplesID)
    finalData$Gene.symbol = data.annotated$Gene.symbol
    ###Mean value if multiple probes
    finalDataWithoutDup = do.call("rbind", splitDF(finalData, nbSplit, samplesID))
    rownames(finalDataWithoutDup) = unique(finalData$Gene.symbol)
    return(finalDataWithoutDup)
}

rmvBatchEffect <- function(dfList, refFile){
    ###Remove batch effect
    set.seed(1234)
    message("Remove batch effect ...")
    data = mergeDf(dfList)
    data = data[, order(names(data))]
    rownames(refFile) = refFile$V3
    keepCohorte = subset(refFile, rownames(refFile) %in% colnames(data))
    excludeCohorte = subset(refFile, !(rownames(refFile) %in% colnames(data)))
    if (dim(excludeCohorte)[1] > 0){
        cohortesExclude = unique(excludeCohorte$V2)
        for (cohorte in cohortesExclude){
            message(paste("No common genes with other cohorte for cohorte : ", cohorte, sep=""))
            message(paste(cohorte, "will be exclude from the analysis", sep=" "))
        }
    }
    keepCohorte = keepCohorte[order(keepCohorte$V3),]
    scaledData = removeBatchEffect(data, batch = keepCohorte$V2)
    return(scaledData)
}

clustering <- function(data, genes, outDir){
    bursteinGenes = c('DHRS2', 'GABRP', 'AGR2', 'PIP', 'FOXA1', 'PROM1', 'NAT1', 'BCL11A', 'ESR1', 'FOXC1', 'CA12', 'TFF3', 'TFF1', 'SCUBE2', 'SFRP1', 'ERBB4','SIDT1', 'PSAT1', 'CHI3L1', 'AR', 'CD36', 'OGN', 'ABCA2', 'CFD', 'IGF1', 'HBB', 'CDH1', 'MEOX2', 'GPX3', 'SCARA5', 'PDK4', 'ENPP2', 'AGTR1', 'LEP', 'LPL', 'DPT', 'TIMP4', 'FHL1', 'SRPX', 'EDNRB', 'SERPINB5', 'SOX10', 'IRX1', 'MIA', 'DSC2', 'TTYH1', 'COL9A3', 'FGL2', 'RARRES3', 'PDE9A', 'BST2', 'PTGER4', 'KCNK5', 'PSMB9', 'HLA-DMA', 'EPHB3', 'IGSF6', 'ST3GAL6', 'RHOH', 'SGPP1','CXCL9', 'CXCL11', 'GBP5', 'GZMB', 'LAMP3', 'GBP1', 'ADAMDEC1', 'CCL5', 'SPON1', 'PBK', 'STAT1', 'EZH2', 'PLAT', 'TAP2', 'SLAMF7', 'HERC5', 'SPOCK1', 'TAP1', 'CD2', 'AIM2')
    dataBurstein = subset(data, rownames(data) %in% bursteinGenes)
    message("Kmeans computing ...")
    res.km = findBestKmeans(dataBurstein, 4, 10000)
    clusters = fviz_cluster(res.km, data = t(dataBurstein),geom = c("point","text"),ellipse.type = "convex", ggtheme = theme_bw()) + ggtitle(paste("Cluster Burstein", dim(dataBurstein)[1], sep=" "))
    dir.create(file.path(outDir, "Clustering"), showWarnings = FALSE)
    QCdir = file.path(outDir, "Clustering")
    ggsave(file.path(outDir, "Clustering", "burstein.jpg"), device="jpeg")
    return(res.km)
}

expressionPlot <- function(data, genes, outDir, res.km){
    conclusion = c()
    dir.create(file.path(outDir, "Boxplot"), showWarnings = FALSE)

    samplesID = colnames(data)

    ovLAR = subset(genes, Subtype == "LAR")$Gene
    ovIM = subset(genes, Subtype == "IM")$Gene
    ovBLIS = subset(genes, Subtype == "BLIS")$Gene
    ovMES = subset(genes, Subtype == "MES")$Gene
    STATS = c("STAT1", "STAT6", "STAT2", "STAT4", "STAT3", "STAT5A", "STAT5B")
    SOX = c("SOX8","SOX18","SOX1","SOX15","SOX5","SOX21","SOX11","SOX14","SOX6","SOX7","SOX1","SOX2","SOX10","SOX12","SOX2","SOX3","SOX4","SOX9","SOX17","SOX13","SOX30")
    osteoAdipo = c("OGN", "ADIPOQ", "PLIN1", "IGF1")
    estrogenRegulated = c("ESR1", "PGR", "FOXA1", "FOXA2", "FOXA3", "GATA3")

    larPlot = meanBoxplot(data, ovLAR, "LAR", samplesID, res.km, conclusion)
    mesPlot = meanBoxplot(data, ovMES, "MES", samplesID, res.km, conclusion)
    blisPlot = meanBoxplot(data, ovBLIS, "BLIS", samplesID, res.km, conclusion)
    imPlot = meanBoxplot(data, ovIM, "IM", samplesID, res.km, conclusion)

    if (!((is.null(larPlot$genePlot)) & (is.null(mesPlot$genePlot)) & (is.null(blisPlot$genePlot)) & (is.null(imPlot$genePlot)))){
        ggarrange(larPlot$genePlot, mesPlot$genePlot, blisPlot$genePlot, imPlot$genePlot, ncol=2, nrow=2) %>% ggexport(filename = file.path(outDir, "Boxplot", "overExpressed.jpg"), width=800, height=900, res = 120)
    }

    IL6plot = meanBoxplot(data, "IL6", "IL6", samplesID, res.km, conclusion)
    JAK1plot = meanBoxplot(data, "JAK1", "JAK1", samplesID, res.km, conclusion)
    STAT3plot = meanBoxplot(data, "STAT3", "STAT3", samplesID, res.km, conclusion)
    if (!((is.null(IL6plot$genePlot)) & (is.null(JAK1plot$genePlot)) & (is.null(STAT3plot$genePlot)))){
        ggarrange(IL6plot$genePlot, JAK1plot$genePlot, STAT3plot$genePlot, ncol=3) %>% ggexport(filename = file.path(outDir, "Boxplot", "IL6_JAK1_STAT3.jpg"), width=1250, height=660, res = 120)
    }

    statsPlot = meanBoxplot(data, STATS, "STAT genes", samplesID, res.km, conclusion)
    if (!(is.null(statsPlot$genePlot))){
        ggarrange(statsPlot$genePlot, ncol=1) %>% ggexport(filename = file.path(outDir, "Boxplot", "STATS.jpg"), width=668, height=662, res = 120)
    }

    soxPlot = meanBoxplot(data, SOX, "SOX genes", samplesID, res.km, conclusion)
    if (!(is.null(soxPlot$genePlot))){
        ggarrange(soxPlot$genePlot, ncol=1) %>% ggexport(filename = file.path(outDir, "Boxplot", "SOX.jpg"), width=668, height=800, res = 120)
    }

    adipoPlot = meanBoxplot(data, osteoAdipo, "Osteocytes and adipocytes specific genes", samplesID, res.km, conclusion)
    if (!(is.null(adipoPlot$genePlot))){
        ggarrange(adipoPlot$genePlot, ncol=1) %>% ggexport(filename = file.path(outDir, "Boxplot", "MES_specific.jpg"), width=668, height=662, res = 120)
    }

    estrogenPlot = meanBoxplot(data, estrogenRegulated, "Estrogen regulated genes", samplesID, res.km, conclusion)
    if (!(is.null(estrogenPlot$genePlot))){
        ggarrange(estrogenPlot$genePlot, ncol=1) %>% ggexport(filename = file.path(outDir, "Boxplot", "LAR_specific.jpg"), width=668, height=662, res = 120)
    }

    FOXC1plot = meanBoxplot(data, c("FOXC1"), "FOXC1", samplesID, res.km, conclusion)
    if (!(is.null(FOXC1plot$genePlot))){
        ggarrange(FOXC1plot$genePlot, ncol=1) %>% ggexport(filename = file.path(outDir, "Boxplot", "FOXC1.jpg"), width=668, height=662, res = 120)
    }

    ARplot = meanBoxplot(data, c("AR"), "AR", samplesID, res.km, conclusion)
    if (!(is.null(ARplot$genePlot))){
        ggarrange(ARplot$genePlot, ncol=1) %>% ggexport(filename = file.path(outDir, "Boxplot", "AR.jpg"), width=668, height=662, res = 120)
    }

    MKIplot = meanBoxplot(data, c("MIB1"), "MKI67", samplesID, res.km, conclusion)
    if (!(is.null(MKIplot$genePlot))){
        ggarrange(MKIplot$genePlot, ncol=1) %>% ggexport(filename = file.path(outDir, "Boxplot", "MKI67.jpg"), width=668, height=662, res = 120)
    }

    data = as.data.frame(t(data))
    data$cluster = res.km$cluster
    write.csv(data, file.path(outDir, "final_clusters.csv"))

    conclusions = as.data.frame(c(imPlot$conclusion, larPlot$conclusion, blisPlot$conclusion, mesPlot$conclusion))
    colnames(conclusions) = "conclusion"
    write.csv(conclusions, file.path(outDir, "conclusion.csv"))
}

formateData <- function(genesListFolder, df){
    types = c("LAR", "IM", "BLIS", "MES")
    finalSubtype = data.frame()
    for (type in types){
        genes = read.table(paste(genesListFolder, type, ".txt", sep = ""), header = FALSE)
        genes$type = rep(type, dim(genes)[1])
        finalSubtype = rbind(finalSubtype, genes)
    }
    colnames(finalSubtype) = c("Gene", "Subtype")

    df = subset(df, rownames(df) %in% finalSubtype$Gene)

    finalExpression = data.frame()
    for (value in rownames(df)){
        tmp = as.data.frame(t(df[value, ]))
        tmp$Gene = rep(value, dim(tmp)[1])
        tmp$Sample = rownames(tmp)
        colnames(tmp) = c("Expression", "Gene", "Sample")
        finalExpression = rbind(finalExpression, tmp)
    }

    return(merge(finalSubtype, finalExpression, by = c("Gene"))[, -1])
}

clusteringGenes <- function(genesListFolder){
    types = c("LAR", "IM", "BLIS", "MES")
    finalSubtype = data.frame()
    for (type in types){
        genes = read.table(paste(genesListFolder, type, ".txt", sep = ""), header = FALSE)
        genes$type = rep(type, dim(genes)[1])
        finalSubtype = rbind(finalSubtype, genes)
    }
    colnames(finalSubtype) = c("Gene", "Subtype")
    return(finalSubtype)
}

main <- function(inputFile, outDir, mode, nbSplit, genesListFolder){
    genes = clusteringGenes(genesListFolder)
    if (mode == "microarray"){
        h = hash()
        cleanedData = list()
        data = read.table(inputFile, sep="\t")
        cohortes = unique(data$V2)
        for (cohorte in cohortes){
            cohorteData = subset(data, V2 == cohorte)
            chipTypes = unique(cohorteData$V4)
            for (chipType in chipTypes){
                chipData = subset(cohorteData, V4 == chipType)
                h[[paste(cohorte, chipType, sep="|")]] = list(data=getData(chipData$V1, chipData$V3), files=chipData$V1)
            }
        }

        for (key in keys(h)){
            chipType = unlist(strsplit(key, "\\|"))[2]
            computeQC(h[[key]]$data, outDir, key, h[[key]]$files)
            cleanedData = append(cleanedData, list(preProcessing(h[[key]]$data, chipType, nbSplit)))
        }

        if (length(unique(data$V2)) > 1 ){
            scaledData = rmvBatchEffect(cleanedData, data)
        }else{
            scaledData = cleanedData[[1]]
        }
    }else{
        scaledData = read.table(inputFile, sep="\t", row.names=1, header=TRUE, check.names = FALSE)
    }
    res.km = clustering(scaledData, genes$Gene, outDir)
    expressionPlot(scaledData, genes, outDir, res.km)
}

option_list = list(
make_option(c("-g", "--geneFilesFolder"), type="character", help="Folder containing genes list files"),
make_option(c("-i", "--inFile"), type="character", help="File containing cel path, cohorte name, sample ID and chip type"),
make_option(c("-m", "--mode"), type="character", help="Data type to cluster. Must be rnaseq or microarray"),
make_option(c("-o", "--outDir"), type="character", help="Output folder"),
make_option(c("-s", "--split"), type="double", help="Number of data split to perform for crossing gene name with probes name. Depending on specs machine, do not perform split can lead to out of memory exception", default=200));
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

main(opt$inFile, opt$outDir, opt$mode, opt$split, opt$geneFilesFolder)
