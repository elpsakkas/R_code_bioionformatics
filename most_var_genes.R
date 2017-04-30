library(DESeq2)
library(statmod)
library(genefilter)


countdata <- read.csv(file.choose(),row.names=1)

countdata <- countdata[rowSums(countdata>0)>5, ]
countdata <- countdata[rowSums(countdata)>1, ]

geneTypes <- factor( c( ENSM="ENSM", ERCC="ERCC" )[
  substr( rownames(countdata), 1, 4 ) ] )


countsMmus <- countdata[ which( geneTypes=="ENSM" ), ]
countsERCC <- countdata[ which( geneTypes=="ERCC" ), ]



sfERCC <- estimateSizeFactorsForMatrix( countsERCC )
sfMmus <- estimateSizeFactorsForMatrix( countsMmus )


nCountsERCC <- t( t(countsERCC) / sfERCC )
nCountsMmus <- t( t(countsMmus) / sfMmus )

meansMmus <- rowMeans( nCountsMmus )
varsMmus <- rowVars( nCountsMmus )
cv2Mmus <- varsMmus / meansMmus^2

meansERCC <- rowMeans( nCountsERCC )
varsERCC <- rowVars( nCountsERCC )
cv2ERCC <- varsERCC / meansERCC^2



par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9)
smoothScatter(log(meansERCC),log(cv2ERCC))

minMeanForFit <- unname( quantile( meansERCC[ which( cv2ERCC > .3 ) ], .95 ) )
useForFit <- meansERCC >= minMeanForFit 
fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansERCC[useForFit] ),cv2ERCC[useForFit] )
a0 <- unname( fit$coefficients["a0"] )
a1 <- unname( fit$coefficients["a1tilde"])
xg <- exp(seq( min(log(meansERCC[meansERCC>0])), max(log(meansERCC)), length.out=1000 ))
vfit <- a1/xg + a0

lines( log(xg), log(vfit), col="black", lwd=3 )
df <- ncol(nCountsERCC) - 1
lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="red")
lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="red")

afit <- a1/meansMmus+a0
varFitRatio <- varsMmus/(afit*meansMmus^2)
varorder <- order(varFitRatio,decreasing=T)


points(log(meansMmus[varorder[1:500]]),log(cv2Mmus[varorder[1:500]]),col='red',cex=2)

legend('bottomleft','Genes (500 top high var)',pch=1,col=c('red'),cex=0.7,bty="n")

pval <- pchisq(varFitRatio*df,df=df,lower.tail=F)
adj.pval <- p.adjust(pval,"fdr")
sig <- adj.pval<1e-3;
table(sig)

log2RelExprMus <- log2(nCountsMmus+1)
highVarTable <- data.frame(
  row.names = NULL,
  ensembl_gene_id = rownames(countsMmus)[ sig ],
  meanNormCount = meansMmus[ sig ],
  # strongest = factor( colnames( log2RelExprMus )[
  #  apply( log2RelExprMus[ sig, ], 1, which.max ) ] ),
  log2RelExprMus[ sig, ],
  check.names=FALSE )

head( highVarTable ) 

write.csv( highVarTable, file = "highly_variant_genes_ALL_filter.csv", row.names=FALSE )
