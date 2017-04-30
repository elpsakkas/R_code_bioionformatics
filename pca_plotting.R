# Created by Asa Björklund 150324
# functions for plotting and running PCA on rpkm-data

run.pca<-function(data,seln=0,log.vals=TRUE,samples.col=TRUE,center=TRUE){
	if (!samples.col){ data<-t(data) }

	# remove zero read genes
	z<-which(rowSums(data)==0)
	if (length(z>0)){
	    data<-data[-z,]
	}

	if (log.vals) { data <- log2(data+0.1) }
	
	# select top varied genes
	if (seln>0){
	    sD<-apply(data,1,sd)
	    o<-order(sD,decreasing=T)
	    data<-data[o[1:seln],]
	}


        myPca <- prcomp(t(data),center=center,scale.=FALSE)
        vars<- myPca$sdev^2
        vars<- vars/sum(vars)
        pcnames<-mat.or.vec(1,length(vars))
        for (i in 1:length(vars)) {
            pcnames[i]<-sprintf("PC%d\n%.5f",i,vars[i])
        }

	myPca$pc.names<-as.vector(pcnames)
        return(myPca)
}
		

 
pca.plot <- function (data, seln=0, selpc=c(1,2),cex=0.7, log.vals=TRUE,center=TRUE,samples.col=TRUE, ...){
	 # data is either a prcomp object or a data matrix

	 if (class(data) != "prcomp"){
	    data<-run.pca(data,seln=seln,log.vals=log.vals,samples.col=samples.col,center=center)
	 }
	 
	 tmpPca <- as.data.frame(data$x[,selpc])
	 colnames(tmpPca)<-data$pc.names[selpc]
    	 plot(tmpPca, ...)
         invisible(data)
}



pca.loadings.plot <- function(data,main="pca loadings",seln=0,log.vals=TRUE,samples.col=TRUE,nPlot=20,nPC=5,center=TRUE,horizontal=TRUE){
	 # data is either a prcomp object or a data matrix

	 if (class(data) != "prcomp"){
	    data<-run.pca(data,seln=seln,log.vals=log.vals,samples.col=samples.col,center=center)
	 }
	 

         aload <- abs(data$rotation[,1:nPC])
         contr<- sweep(aload, 2, colSums(aload), "/")
         par(mfrow=c(nPC,2),mar=c(1,10,2,1),oma=c(1,5,1,1),cex=0.4)
         for (i in 1:nPC) {
             top<-order(data$rotation[,i],decreasing=T)[1:nPlot]
             bottom<-order(data$rotation[,i],decreasing=F)[1:nPlot]
             barplot(contr[top,i],main=sprintf("genes on pos axis PC%d",i),ylab="% contr",las=2,horiz=horizontal)
             barplot(contr[bottom,i],main=sprintf("genes on neg axis PC%d",i),ylab="% contr",las=2,horiz=horizontal)
         }
         invisible(data)
}
    
