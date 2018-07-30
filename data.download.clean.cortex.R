########## cader group brainspan wgcna ########

setwd('....')
download.file('http://www.brainspan.org/api/v2/well_known_file_download/267666525', destfile = 'brainspan.r.zip')
unzip('brainspan.r.zip',exdir='brainspan')

rwz <- read.csv('~/brainspan/rows_metadata.csv', stringsAsFactors = F)
cls <- read.csv('~/brainspan/columns_metadata.csv', stringsAsFactors = F)
exprs <- read.csv('~/brainspan/expression_matrix.csv', row.names = 1, header = F)

#### only cortical samples #####

cls.c <- cls[grep('cortex', cls$structure_name, ignore.case = T), ]

exprs.c <- exprs[, cls.c$column_num]

tc <- exprs

tc <- ifelse(exprs.c == 0, T, F)

zeros <- rowSums(tc)
# 391 * .8 = 312.8
# 391 - 313  -> 78
# 170 *. 8  --> 170 - 136 34

zeros.filter <- zeros <= 34
# filters for where zeros are more  than 20% of samples

exprs.c <- exprs.c[zeros.filter, ]
rwz.c <- rwz[zeros.filter, ]

exprs.c <- exprs.c[, grep('pcw', cls.c$age)]
cls.c <- cls.c[grep('pcw', cls.c$age), ]


exprs.log <- log2(exprs.c + 1)

### next protein coding????
######

library("PoiClaClu")
poisd <- PoissonDistance(t(exprs.c))

samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste(cls.c$age, cls.c$gender,cls.c$donor_id, sep = '-')
colnames(samplePoisDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

brainspan.exp <- SingleCellExperiment(assays = list(counts = as.matrix(exprs.c), logcounts = as.matrix(exprs.log)))
rwz.c <- DataFrame(rwz.c)
rowData(brainspan.exp) <- rwz.c
cls.c <- DataFrame(cls.c)

rownames(cls.c) <- NULL
colData(brainspan.exp) <- cls.c

scater::plotPCA(brainspan.exp, ncomponents = 4, colour_by = 'structure_name')

