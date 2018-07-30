cortecon.txi <- readRDS('~/Desktop/allen/brainspan/cortecon/txi.cort.alone.RDS')
sampletable <- readRDS('~/Desktop/allen/brainspan/cortecon/sample.info.cortecon.RDS')
class(sampletable$Characteristics..Stage.)
sampletable$Characteristics..Stage.[sampletable$Characteristics..Stage. == 'Neural Differentation'] <- 'Neural Differentiation'
sampletable$Characteristics..Stage. <- factor(sampletable$Characteristics..Stage., levels = c('Pluripotency', 
                                                                                              'Neural Differentiation',
                                                                                              'Cortical Specification',
                                                                                              'Deep Layer Formation',
                                                                                              'Upper Layer Formation'))

cortecon.data <- DESeqDataSetFromTximport(cortecon.txi, sampletable, design = ~Characteristics..Stage.)  

cortecon.data <- cortecon.data[rowSums(counts(cortecon.data)) > 1, ]

vsd <- vst(cortecon.data)
rld <- rlog(cortecon.data)

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$Characteristics..Stage., vsd$Comment..Sample_title., vsd$Scan.Name, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)#, file = 'rlog.pheatmap.png')

#####
poisd <- PoissonDistance(t(counts(cortecon.data)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( vsd$Characteristics..Stage., vsd$Comment..Sample_title., sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

plotPCA(rld, intgroup = c('samp'))

test <- SingleCellExperiment(assays = list(counts = counts(cortecon.data), logcounts = assay(rld)))
rowData(test) <- rowData(rld)
colData(test) <- colData(rld)
scater::plotPCA(test, ncomponents = 4, colour_by = 'samp')

dds <- DESeq(cortecon.data)


#######################
rld <- readRDS('cortecon/cortecon.normalized.data.rlog.RDS')
wgcna with rld 
wgcna.cort <- assay(rld)

sampleTree = hclust(dist(wgcna.cort), method = 'average')
saveRDS(rld, file = 'cortecon.normalized.data.rlog.RDS')
saveRDS(dds, file = 'cortecon.deg.RDS')

bcftools index -t $n
