library("DESeq2")

countdata <- read.csv(file.choose(),row.names=1)

coldata <- read.csv(file.choose(),row.names=1)

countdata <- countdata[rowSums(countdata>0)>1, ]
countdata <- countdata[rowSums(countdata)>1, ]

geneTypes <- factor( c( ENSM="ENSM", ERCC="ERCC" )[
  substr( rownames(countdata), 1, 4 ) ] )

countsMmus <- countdata[ which( geneTypes=="ENSM" ), ]
countsERCC <- countdata[ which( geneTypes=="ERCC" ), ]

# sfERCC <- estimateSizeFactorsForMatrix( countsERCC )
sfMmus <- estimateSizeFactorsForMatrix( countsMmus )

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~ clusterboot)

# sizeFactors(dds) <- sfERCC
sizeFactors(dds) <- sfMmus


vsd <- varianceStabilizingTransformation(dds)

plotPCA(vsd, intgroup=c("clusterboot"))


dds <- DESeq(dds)
res <- results(dds, alpha=0.05)

sum(res$padj < 0.05, na.rm=TRUE)

resSig <- subset(res, padj < 0.05)

write.csv(as.data.frame(resSig),file="DE_club_cilia.csv")


##########################################################
# Optimize the plotPCA function
library(genefilter)
ntop <- 500
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- ( assay(vsd)[select, ] )
## Continue with the PCA function


