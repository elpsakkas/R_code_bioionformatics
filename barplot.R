
clus.2.TF.expr <- rowMeans(imp.TF.all.clusters.merge[,clus_1])



barplot(imp.TF.all.clusters.merge.1$cluster_5_TF,
        width =5,main="Cluster_5_TF",
        xlab="Important genes",
        ylab="LogNormExpression",
        names.arg=imp.TF.all.clusters.merge.1$external_gene_name,
        col="cyan",
        horiz=F,axis.lty=1,cex.names=0.4)
