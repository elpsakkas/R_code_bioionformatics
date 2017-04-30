library("DESeq2")
library("genefilter")
library("statmod")




countdata <- read.csv(file.choose(),row.names=1)
meta <- read.csv(file.choose(),row.names=1)
cols <- meta[,'cell']
cols <- as.character(cols)
countdata <- countdata[,cols]


countdata <- countdata[rowSums(countdata>0)>10, ]
countdata <- countdata[rowSums(countdata)>1, ]

geneTypes <- factor( c( ENSM="ENSM", ERCC="ERCC" )[
  substr( rownames(countdata), 1, 4 ) ] )


countsMmus <- countdata[ which( geneTypes=="ENSM" ), ]
countsERCC <- countdata[ which( geneTypes=="ERCC" ), ]


sfERCC <- estimateSizeFactorsForMatrix( countsERCC )
sfMmus <- estimateSizeFactorsForMatrix( countsMmus )


nCountsERCC <- t( t(countsERCC) / sfERCC )
nCountsMmus <- t( t(countsMmus) / sfMmus )



colERCC <- "#00207040"
colMmus <- "black"


meansMmus <- rowMeans( nCountsMmus )
varsMmus <- rowVars( nCountsMmus )
cv2Mmus <- varsMmus / meansMmus^2

meansERCC <- rowMeans( nCountsERCC )
varsERCC <- rowVars( nCountsERCC )
cv2ERCC <- varsERCC / meansERCC^2

minMeanForFitA <- unname( quantile( meansERCC[ which( cv2ERCC > .3 ) ], .95 ) )
minMeanForFitA
useForFitA <- meansERCC >= minMeanForFitA
fitA <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansERCC[useForFitA] ),
                    cv2ERCC[useForFitA] )

fitA$coefficients


xi <- mean( 1 / sfERCC )
a0 <- unname( fitA$coefficients["a0"] )
a1 <- unname( fitA$coefficients["a1tilde"] - xi )
c( a0, a1 )

plot( NULL, xaxt="n", yaxt="n",
      log="xy", xlim = c( 1e-1, 3e5 ), ylim = c( .005, 8 ),
      xlab = "read counts", ylab = "coefficient of variation (CV^2)" )
axis( 1, 10^(-1:5), c( "0.1", "1", "10", "100", "1000",
                       expression(10^4), expression(10^5) ) )
axis( 2, 10^(-2:1), c( "0.01", "0.1", "1", "10" ), las=2 )
#abline( h=10^(-2:1), v=10^(-1:5), col="#D0D0D0", lwd=2 )

points( meansERCC, cv2ERCC, pch=20, cex=.4, col="blue" )

xg <- 10^seq( -2, 6, length.out=1000 )
lines( xg, (xi+a1)/xg + a0, col="blue", lwd=3 )

df <- ncol(countsMmus) - 1
lines( xg, ( (xi+a1)/xg + a0 ) * qchisq( .975, df ) / df,
       col="#FF000080", lwd=2, lty="dashed" )
lines( xg, ( (xi+a1)/xg + a0 ) * qchisq( .025, df ) / df,
       col="#FF000080", lwd=2, lty="dashed" ) 

psia1theta <- mean( 1 / sfMmus ) + a1 * mean( sfERCC / sfMmus )

minBiolDisp <- .5^2
m <- ncol(countsMmus)
cv2th <- a0 + minBiolDisp + a0 * minBiolDisp
testDenom <- ( meansMmus * psia1theta + meansMmus^2 * cv2th ) / ( 1 + cv2th/m )
p <- 1 - pchisq( varsMmus * (m-1) / testDenom, m-1 )

padj <- p.adjust( p, "BH" )
sig <- padj < .1
sig[is.na(sig)] <- FALSE
table( sig )

plot( NULL, xaxt="n", yaxt="n",
      log="xy", xlim = c( 1e-1, 3e5 ), ylim = c( .005, 8 ),
      main= 'Dispersion plot',
      xlab = "read counts", ylab = "CV^2" )
axis( 1, 10^(-1:5), c( "0.1", "1", "10", "100", "1000",
                       expression(10^4), expression(10^5) ) )
axis( 2, 10^(-2:1), c( "0.01", "0.1", "1", "10" ), las=2 )
#abline( h=10^(-2:1), v=10^(-1:5), col="blue", lwd=2 )

points( meansMmus, cv2Mmus, pch=20, cex=.2,
        col = ifelse( padj < .1, "red", colMmus ) )

xg <- 10^seq( -2, 6, length.out=1000 )
lines( xg, (xi+a1)/xg + a0, col="blue", lwd=2 )

lines( xg, psia1theta/xg + a0 + minBiolDisp, lty="dashed", col="red", lwd=3 )
points( meansERCC, cv2ERCC, pch=20, cex=.4, col="blue" )


# log2RelExprMus <- log2( nCountsMmus / meansMmus )
log2RelExprMus <- log2( nCountsMmus + 1)
highVarTable <- data.frame(
  row.names = NULL,
  ensembl_gene_id = rownames(countsMmus)[ sig ],
  meanNormCount = meansMmus[ sig ],
  #strongest = factor( colnames( log2RelExprMus )[
  #apply( log2RelExprMus[ sig, ], 1, which.max ) ] ),
  log2RelExprMus[ sig, ],
  check.names=FALSE )

head( highVarTable ) 

write.csv( highVarTable, file = "highly_variant_genes_all_417.csv", row.names=FALSE )
