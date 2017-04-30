#########################

#### File with various task functions 
### Author : Akis Sakkas
## Date: 2016

count.genes <- function(x, lim=1){
  length(which(x>lim))
}
                                      
color.data.frame <- function(x,){
  x$Color="black"
  x$Color
  
}
count.genes <- function(x, lim=1){
  length(which(x>lim))
}

count.gfp <- function(x){
  gfp <- x[which(rownames(x)=="pAcGFP1-N1"),]
  t(log2(gfp+1))
}


color.stage <- function(x){
  color.palette <-colorRampPalette(c("cyan1",
                                     "chartreuse1",
                                     "seagreen",
                                     "skyblue3",
                                     "darkblue"),space="Lab")
  color.palette(length(x$Stage))[rank(x$Stage)]
}                   
                   
colStage <- color.stage(met)


legend('topright',c('E16','E19','P2','P21','P60'),pch=c(15,15,15,15,15),
       col=c('cyan1','chartreuse1','seagreen','skyblue3','darkblue'),cex=0.7)
#col.list <- c("cyan1","chartreuse1","seagreen","skyblue3","darkblue")
#palette(col.list)

color.cluster <- function(x){
          coloring <-colorRampPalette(c("red",
                                     "black"),space="Lab")
  coloring(length(number.genes))[rank(number.genes)]
}                   

colCluster <- color.cluster(met)


######### Another way of coloring
  color.palette <- c("red",
                      "orange",
                      "yellow",
                      "green",
                      "cyan",
                      "lightblue",
                      "black")
  
palette(color.palette) # Use col in the function

par=(cex=3.0)
plot(dif,pch = c(2, 0, 20)[as.numeric(met$shape)],
col = c("orange","black","darkgreen","cyan","darkred","grey","blue","red")[as.numeric(meta$clusterboot)]



legend('bottomright',c("cluster_1","cluster_2","cluster_3","cluster_4","cluster_5","cluster_6","cluster_7","cluster_8"),
       col = c("orange","black","darkgreen","cyan","darkred","grey","blue","red"),
       pch=c(20,20,20,20,20,20,20),bty="n",cex=1)
legend('bottomright',c('E16','E19','P2','P21','P60'),pch=c(2,0,0,20,20),,cex=2,bty="n")

## This is the original labeling
plot(tsne_out$Y[,1],tsne_out$Y[,2],pch = c(4,15,16,17,9)[as.numeric(meta$Stage)],xlab='tSNE_1',ylab='tSNE_2')
plot(base.pca$x[,1],base.pca$x[,2],pch = c(4,15,16,17,9)[as.numeric(meta$Stage)],xlab='PC1',ylab='PC2')
legend('bottomright',c('E16','E19','P2','P21','P60'),pch=c(4,15,16,17,9),,cex=1,bty="n")

#################################################





heatmap <- heatmap.2(as.matrix(most.imp.genes),cexRow=0.6,cexCol=0.40,key='True',
                     main= "RF most important genes",
                     scale="none",trace='none',
                     Colv=as.dendrogram(pfit),dendrogram="column",
                     ColSideColor=c("grey","black","darkgreen","darkred","blue","orange","red")[as.numeric(meta$clusterboot)],
                     labCol=F,col=my_palette)

my_palette <- colorRampPalette(c("red", "orange","yellow","white"))(n = 20)

##### Running ggplot2 package for visualization - TSNE

p1 <- ggplot(data=as.data.frame(tsne_out$Y), 
             aes(x=tsne_out$Y[,1],y=tsne_out$Y[,2]))
p1 + geom_point(aes(colour = sox2.raw)) + 
scale_colour_gradient(low = "black",high = 'orange') +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

####Violin plots 1
mat2 <- mat1[,c('Bpfib1',)]
df.m <- reshape2::melt(mat2, id.vars = NULL)
colnames(df.m) <- c('Gene','Expression')
ggplot(df.m, aes(x = Gene, y = Expression)) +  geom_violin(fill = 'grey80', colour = "#3366FF") 
####Violin plots 2 (automated)
plotViolFunc <- function(x, na.rm = TRUE, ...) {
  nm <- colnames(x)
  for (i in seq_along(nm)) {
    print(ggplot(x,aes_string(groups,nm[i])) + geom_violin() + geom_jitter(height = 0)) 
    }
}
pdf('viol_plots.pdf')
plotViolFunc(as.data.frame(mat1))
dev.off()
######################################################    
plotTsneFunc <- function(x, na.rm = TRUE, ...) {
  nm <- list(x)
  for (i in seq_along(nm)) {
    print(ggplot(data=as.data.frame(tsne_out$Y),aes(x=tsne_out$Y[,1],y=tsne_out$Y[,2])) + geom_point(aes(colour=gfp)))
  }
}

plotTsneFunc(mat)

################################################3
#Normalization function DESeq2

calc_sf <-
  function (expr_mat, spikes = NULL) 
  {
    geomeans <- exp(rowMeans(log(expr_mat[-spikes, ])))
    SF <- function(cnts) {
      median((cnts/geomeans)[(is.finite(geomeans) & geomeans > 
                                0)])
    }
    norm_factor <- apply(expr_mat[-spikes, ], 2, SF)
    return(t(t(expr_mat)/norm_factor))
  }
######################################################3

top.loading.5.pc <- function(x, na.rm = TRUE, ...) {
  pcs <- 5
  for (i in 1:pcs) {
    pc[i] <- sort(x$rotation[,i])
    pc[i].head <- as.data.frame(head(pc[i],n=50L))
    pc[i].tail <- as.data.frame(tail(pc[i],n=50L))
    colnames(pc[i].head) <- 'pc'
    colnames(pc[i].tail) <- 'pc'
    pc[i].total <- rbind(pc[i].head,pc[i].tail)
  }
  
}
