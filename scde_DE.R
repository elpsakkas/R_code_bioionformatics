library(scde)
library(samr)


cd <- read.csv(file.choose(),row.names=1)
met <- read.csv(file.choose(),row.names=1)
groups <- met[,'discover']
cols <- met[,'cell']
cols <- as.character(cols)


cd <- cd[,cols]
cd <- cd[rowSums(cd>0)>2, ]
cd <- cd[rowSums(cd)>1, ]
geneTypes <- factor( c( ENSM="ENSM", ERCC="ERCC" )[
   substr( rownames(cd), 1, 4 ) ] )

countsMmus <- cd[ which( geneTypes=="ENSM" ), ]
# sfMmus <- estimateSizeFactorsForMatrix( countsMmus )
# nCountsMmus <- t( t(countsMmus) / sfMmus )

# mat <- nCountsMmus
counts <- countsMmus
# counts<-apply(counts,2,function(x) {storage.mode(x) <- 'integer'; x})

scde.fitted.model <- scde.error.models(counts=counts,groups=groups,save.model.plots=F)
scde.prior <- scde.expression.prior(models=scde.fitted.model,counts=counts)

dif <- scde.expression.difference(scde.fitted.model,counts,scde.prior,groups=groups,
                                  n.randomizations  =  100, n.cores  =  2, verbose  =  1)
###### We can sort the dif object according to Z scores #######

head(dif[order(dif$Z, decreasing  =  TRUE), ]) # top upregulated genes
tail(dif[order(dif$Z, decreasing  =  TRUE), ]) # top downregulated genes

write.csv(dif[order(abs(dif$Z), decreasing = TRUE), ], file = "DE_scde_basc_at2.csv")

scde.test.gene.expression.difference("GENE_NAME", models = o.ifm, counts = cd, prior = o.prior)


####################################################################################
p.values <- 2*pnorm(abs(dif$Z),lower.tail=F) # 2-tailed p-value
p.values.adj <- 2*pnorm(abs(dif$cZ),lower.tail=F) # Adjusted to control for FDR
resSig <- which(p.values.adj<0.05)
length(resSig)



ordered.genes <- order(p.values.adj[resSig]) # order by p-value
de_scde <- cbind(dif[resSig,1:3],p.values.adj[resSig])[ordered.genes,]
colnames(de_scde) <- c("Lower limit","log2 fold change","Upper limit","p-value")

write.csv(de_scde,file="DE_scde_club1_club2.csv")
##############################################################
## SAMseq DE
# counts <- as.matrix(counts)
res <- SAMseq(counts, groups, resp.type = "Two class unpaired",
                 genenames=rownames(counts),fdr.output = 0.05) 


results.table <- rbind(res$siggenes.table$genes.up, res$siggenes.table$genes.lo)
