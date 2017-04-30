library(gplots)


x <-read.table("highly_variant_genes_E19.csv",header=TRUE,row.names=1,sep=",")



x <- x[order(-x$meanNormCount),]

# x <- subset(x, meanNormCount > 550 & meanNormCount < 650)

head(x)

# is.na(x) <- sapply(x, is.infinite)

# x <- na.omit(x)


y <- x[c(0,3:53)]# We choose specific cols that we want to convert them to a matrix

head(y)

y <- y[700:800,]


write.csv(y,file="P60_high_var_top.csv")

w <- read.table("E19_high_var_top.csv",header=TRUE,row.names=1,sep=",")
# w <- read.csv("all_high_var_top.csv",row.names=1)


head(w)

# rownames(w) <- w$geneID # We use the first col as an index column
# is.na(w) <- sapply(w, is.infinite)

# data <- na.omit(w)

w <- as.matrix(w[,-1]) # We scale the values to make them comparable


color.palette = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"), space="Lab")

### And make the heatmap
heatmap.2(w,cexRow=0.45,cexCol=0.40,key=TRUE,density="none",scale="none",trace="none",col="heat.colors",ColSideColor=colcell_type)




