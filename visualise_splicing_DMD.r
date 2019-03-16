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
sample_name <- paste(dirname, "dmd", sep="/")

si_dmd = DataFrame(sample_name = sample_name, file_bam = c("dmd.sorted.bam"), paired_end = FALSE, read_length=96, frag_length=200, lib_size=10978)
si_dmd

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
txdb <- keepSeqlevels(txdb, "chrX")
seqlevelsStyle(txdb) <- "UCSC"

gr_dmd <- GRanges(seqnames="chrX", ranges=IRanges(start=31097677, end=33339441), strand="-")

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

