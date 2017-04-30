library(plotly)
library(DESeq2)

cd <- read.csv(file.choose(),row.names=1)
met <- read.csv(file.choose(),row.names=1)

cd <- cd[rowSums(cd>0)>10, ]
cd <- cd[rowSums(cd)>1, ]

geneTypes <- factor( c( ENSM="ENSM", ERCC="ERCC" )[
  substr( rownames(cd), 1, 4 ) ] )

countsMmus <- cd[ which( geneTypes=="ENSM" ), ]
countsERCC <- cd[ which( geneTypes=="ERCC" ), ]

sfERCC <- estimateSizeFactorsForMatrix( countsERCC )
# sfMmus <- sfERCC
sfMmus <- estimateSizeFactorsForMatrix( countsMmus )


nCountsERCC <- t( t(countsERCC) / sfERCC )
nCountsMmus <- t( t(countsMmus) / sfMmus )

LogNcountsMmus=log2(nCountsMmus+1)

# met$Stage <- factor(met$Type,labels = c("single-read", "paired-end"))




number.genes <- apply(LogNcountsMmus,2,count.genes)



plot_ly(base.pca, 
        x = base.pca$x[,1],
        y = base.pca$x[,2],
        text= meta$cell, 
        mode = "markers")
        #color = met_pseudo$Stage)
        #size = ratio)

##############################
library(FactoMineR)
library(factoextra)

fviz_pca_biplot(base.pca, label="var",
                #habillage=met$hc_cluster,
                addEllipses=TRUE,
                ellipse.level=0.95,
                axes = c(1,3),
                select.var = list(contrib=10,labelsize=0.3))

##############################
library(Rtsne)
d <- stats::dist(t(mat))
set.seed(123) 
tsne_out <- Rtsne(d,
                  dims=2,
                  initial_dims=5,
                  is_distance=TRUE, 
                  perplexity=32, 
                  verbose = 2,
                  theta=0.1,pca=T) 
####tsne - 2d
plot_ly(#tsne_out,
        x=tsne_out$Y[,1],
        y =tsne_out$Y[,2],
        mode="markers", 
        text=meta$cell, 
        color=spc, type='scatter')


#####tsne - 3d
plot_ly(
        x=tsne_out$Y[,1],
        y =tsne_out$Y[,2],
        z =tsne_out$Y[,3],
        type="scatter3d",
        mode="markers", 
        text=meta$cell, 
        color=meta$Stage)


###############################
library(fastICA)

plot_ly(ics, 
        x = ics$A[1,],
        y = ics$A[2,],
        text= met$Run, 
        mode = "markers",
        color = met$Stage,
        size = log.GFP[,1])

##############################
library(destiny)

plot_ly(dif,
        x = dif$eigenvectors[,1],
        y = dif$eigenvectors[,2],
        text= met$Run, 
        mode = "markers",
        color = met$Stage)

