library(MAST)
library(DESeq2)
library(data.table)

cd <- read.csv(file.choose(),row.names=1)
met <- read.csv(file.choose(),row.names=1,sep=row.names)
#met <- met[met$clusterboot %in% c('cluster_5','cluster_7'),]

cols <- met[,'cell']
cols <- as.character(cols)

cd <- cd[,cols]
cd <- cd[rowSums(cd>0)>2, ]
cd <- cd[rowSums(cd)>1, ]

met$rf <- 'unknown'
met$rf <- ifelse(met$clusterboot=='cluster_5',met$rf,2)
groups <- as.factor(met[,'rf'])



geneTypes <- factor( c( ENSM="ENSM", ERCC="ERCC" )[
  substr( rownames(cd), 1, 4 ) ] )

countsMmus <- cd[ which( geneTypes=="ENSM" ), ]

sfMmus <- estimateSizeFactorsForMatrix( countsMmus )

norm.CountsMmus <- t( t(countsMmus) / sfMmus )

Log.counts=log2(norm.CountsMmus+1) + 1

mat <- Log.counts

feat <- rownames(mat)

feat <- as.data.frame(feat)

colnames(feat) <- c("primerid")



# write.csv(feat, file= 'fdata.csv')
# feat <- read.csv('fdata.csv',row.names=1)

scmat <- FromMatrix(mat, met, feat) 

zlm.groups <- zlm.SingleCellAssay(~groups, scmat) # create the model

summary.groups <- summary(zlm.groups)
summary.groups

summary.groups <- summary(zlm.groups, doLRT='groupsunknown') 

summary.all <- summary.groups$datatable

fc.table <- merge(summary.all[contrast=='groupsunknown' & component=='H',
                              .(primerid, `Pr(>Chisq)`)], 
                  summary.all[contrast=='groupsunknown' & component=='logFC',
                              .(primerid, coef, ci.hi, ci.lo)], by='primerid')

fc.table[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]

#fc.table.Sig <- merge(fc.table[fdr<.05 & abs(coef)>0],
#                     as.data.table(mcols(scmat)), by='primerid')

fc.table.Sig <- merge(fc.table,
                      as.data.table(mcols(scmat)), by='primerid')

colnames(fc.table.Sig) <- c('ensembl_gene_id','Pr>Chisq','LogFC','ci.hi','ci.lo','fdr')

write.csv(fc.table.Sig, file = 'MAST_DE_basc1_vs_basc2.csv')


