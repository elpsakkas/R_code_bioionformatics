source("functions.R")
library(DESeq2)
library(gplots)


cd <- read.csv("counts_ALL_GFP.csv",row.names=1)
met <- read.csv("metadata_all_samples.csv",row.names=1)

cd <- cd[rowSums(cd>0)>10, ]
cd <- cd[rowSums(cd)>10, ]

geneTypes <- factor( c( ENSM="ENSM", ERCC="ERCC" )[
  substr( rownames(cd), 1, 4 ) ] )

countsMmus <- cd[ which( geneTypes=="ENSM" ), ]
countsERCC <- cd[ which( geneTypes=="ERCC" ), ]

sfERCC <- estimateSizeFactorsForMatrix( countsERCC )
sfMmus <- sfERCC

nCountsERCC <- t( t(countsERCC) / sfERCC )
nCountsMmus <- t( t(countsMmus) / sfMmus )

Log.counts.Mmus=log2(nCountsMmus+1)


C <- apply(Log.counts.Mmus,2,cor,Log.counts.Mmus)
# diag(C)<-NA
rownames(C) <- colnames(Log.counts.Mmus)

color.palette = colorRampPalette(c("midnightblue","dodgerblue3",
                                 "white","goldenrod1","darkorange2"),
                                 space="Lab")


colBatch <- make_colors(met$Batch,unique(met$Batch),rainbow(3))
# colStage <- make_colors(met$Stage,unique(met$Stage),rainbow(5))

### And make the heatmap
heatmap <- heatmap.2(mat_x,cexRow=0.35,cexCol=0.35,key=TRUE,
           density="none",main= "",
           scale="none",
           trace="none")
           col=color.palette)
           ColSideColor=pfit)
#           RowSideColor=colStage)

####################### Coloring of the developmental stages
color.palette = colorRampPalette(c("aliceblue",
                                   "khaki2",
                                   "aquamarine3",
                                   "cyan4"),space="Lab")
 #### Maybe?
colStage <- color.palette(length(met$Stage))[rank(met$Stage)]
