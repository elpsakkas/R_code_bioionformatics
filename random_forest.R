###### Random forest permutations with importance features



library(randomForest)
library(rfPermute)
library(DESeq2)

start.time <- Sys.time()



#cd <- read.csv('/home/rstudio/SingleCellProject/R_wdir/Final_data/counts_filtered_final_441.csv',row.names=1)
#meta <- read.csv('/home/rstudio/SingleCellProject/R_wdir/Final_data/metadata_all_filtered_final_441.csv',row.names=1)
#load('matrices_analysis_441.RData')





#cd <- cd[rowSums(cd>0)>10, ]
#cd <- cd[rowSums(cd)>10, ]

#geneTypes <- factor( c( ENSM="ENSM", ERCC="ERCC" )[
#  substr( rownames(cd), 1, 4 ) ] )


#countsERCC <- cd[ which( geneTypes=="ERCC" ), ]
#sfERCC <- estimateSizeFactorsForMatrix( countsERCC )

#countsMmus <- cd[ which( geneTypes=="ENSM" ), ]
#sfMmus <- estimateSizeFactorsForMatrix( countsMmus )
#nCountsMmus <- t( t(countsMmus) / sfMmus )
#LogNcountsMmus=log2(nCountsMmus+1)

#mat <- LogNcountsMmus
meta$rf <- 'unknown'
meta$rf <- ifelse(meta$clusterboot2=='1',meta$rf,2)
groups <- meta$rf


cluster_1 <- rfPermute(x = t(mat), y = as.factor(groups), ntree = 1000, 
                        importance = T, proximity = T, 
                        do.trace = T, varImpPlot = T,
                        varUsed = T,nrep = 0, num.cores = 1)

meta$rf <- 'unknown'
meta$rf <- ifelse(meta$clusterboot2=='2',meta$rf,2)
groups <- meta$rf


cluster_2 <- rfPermute(x = t(mat), y = as.factor(groups), ntree = 1000, 
                       importance = T, proximity = T, 
                       do.trace = T, varImpPlot = T,
                       varUsed = T,nrep = 0, num.cores = 1)


meta$rf <- 'unknown'
meta$rf <- ifelse(meta$clusterboot2=='3',meta$rf,2)
groups <- meta$rf


cluster_3 <- rfPermute(x = t(mat), y = as.factor(groups), ntree = 1000, 
                       importance = T, proximity = T, 
                       do.trace = T, varImpPlot = T,
                       varUsed = T,nrep = 0, num.cores = 1)

meta$rf <- 'unknown'
meta$rf <- ifelse(meta$clusterboot2=='4',meta$rf,2)
groups <- meta$rf


cluster_4 <- rfPermute(x = t(mat), y = as.factor(groups), ntree = 1000, 
                       importance = T, proximity = T, 
                       do.trace = T, varImpPlot = T,
                       varUsed = T,nrep = 0, num.cores = 1)



meta$rf <- 'unknown'
meta$rf <- ifelse(meta$clusterboot2=='5',meta$rf,2)
groups <- meta$rf


cluster_5 <- rfPermute(x = t(mat), y = as.factor(groups), ntree = 1000, 
                       importance = T, proximity = T, 
                       do.trace = T, varImpPlot = T,
                       varUsed = T,nrep = 0, num.cores = 1)


meta$rf <- 'unknown'
meta$rf <- ifelse(meta$clusterboot2=='6',meta$rf,2)
groups <- meta$rf


cluster_6 <- rfPermute(x = t(mat), y = as.factor(groups), ntree = 1000, 
                       importance = T, proximity = T, 
                       do.trace = T, varImpPlot = T,
                       varUsed = T,nrep = 0, num.cores = 1)

meta$rf <- 'unknown'
meta$rf <- ifelse(meta$clusterboot2=='7',meta$rf,2)
groups <- meta$rf


cluster_7 <- rfPermute(x = t(mat), y = as.factor(groups), ntree = 1000, 
                       importance = T, proximity = T, 
                       do.trace = T, varImpPlot = T,
                       varUsed = T,nrep = 0, num.cores = 1)








#plot(rp.importance(cluster_1), alpha = 0.05, sig.only = T,n=50,type="MeanDecreaseGini")

#write.csv(rp.importance(cluster_1),file="important_genes_cluster_1.csv")

p1 <- impHeatmap(cluster_1,n=25,ranks = F)
p2 <- impHeatmap(cluster_2,n=25,ranks = F)
p3 <- impHeatmap(cluster_3,n=25,ranks = F)
p4 <- impHeatmap(cluster_4,n=25,ranks = F)
p5 <- impHeatmap(cluster_5,n=25,ranks = F)
p6 <- impHeatmap(cluster_6,n=25,ranks = F)
p7 <- impHeatmap(cluster_7,n=25,ranks = F)




end.time <- Sys.time()
time.taken <- end.time - start.time

time.taken

#save.image('RF_clusters.RData')
