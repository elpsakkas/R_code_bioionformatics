library("Rsamtools")
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")


dir <- "/Volumes/My Passport for Mac/Treutlein_E18/"

sampleTable <- read.csv(file.choose(),row.names=1)

filenames <- file.path(dir, paste0(sampleTable$cell, ".merge.sort.bam"))

file.exists(filenames)
bamfiles <- BamFileList(filenames, yieldSize=2000000)

gtffile <- file.path(dir,"ref-transcripts-pAcGFP1-N1-ERCC_spikes.gtf")

txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character())
ebg <- exonsBy(txdb, by="gene")

se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=T,
                        ignore.strand=TRUE)

counts <- assay(se)
write.csv(counts, file="counts_Treutlein_E18.csv")





