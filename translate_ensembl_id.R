library("biomaRt")

listMarts(host="www.ensembl.org")



ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                  dataset="mmusculus_gene_ensembl", 
                  host = "oct2016.archive.ensembl.org")


cd <- read.csv(file.choose(), header =T)
head(cd)
genes <- cd$ensembl_gene_id

gene.list <- getBM(filters= "ensembl_gene_id",
                             attributes= c("ensembl_gene_id",
                                           "external_gene_name"),
                                           values=genes,
                                           mart= ensembl)

gene.name <- merge(cd,gene.list, by = "ensembl_gene_id")

write.csv(gene.name, file = "highly_variant_genes_all_442.csv")



