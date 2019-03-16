# Install the package
source("https://bioconductor.org/biocLite.R")
biocLite("SGSeq")
biocLite("TxDb.Hsapiens.UCSC.hg38.knownGene")

# Load data
library(SGSeq)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

dirname <- "/home/picrin/programming/muscle_splicing"
setwd(dirname)
getwd()
list.files()
sample_name <- paste(dirname, "serpine1", sep="/")

si_dmd = DataFrame(sample_name = sample_name, file_bam = c("serpine1.sorted.bam"), paired_end = FALSE, read_length=96, frag_length=200, lib_size=23109)
si_dmd

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
txdb <- keepSeqlevels(txdb, "chr7")
seqlevelsStyle(txdb) <- "UCSC"
  
gr_dmd <- GRanges(seqnames="chr7", ranges=IRanges(start=101127089, end=101139266), strand="+")

txf_ucsc <- convertToTxFeatures(txdb)
txf_ucsc <- txf_ucsc[txf_ucsc %over% gr_dmd]
head(txf_ucsc)

sgfc_ucsc <- analyzeFeatures(si_dmd, features = txf_ucsc)
sgfc_ucsc

type(txf_ucsc)

head(txName(txf_ucsc))

colData(sgfc_ucsc)
rowRanges(sgfc_ucsc)
head(counts(sgfc_ucsc))
head(FPKM(sgfc_ucsc))

df <- plotFeatures(sgfc_ucsc, geneID = 1)

