library("knitr")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BiocStyle")
library("BiocStyle")
.cran_packages <- c("ggplot2", "gridExtra")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(.bioc_packages[!.inst], ask = F)
}
if(any(!.inst)) {
  BiocManager::install(.bioc_packages)
}
BiocManager::install("dada2", version = "3.9")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")












library(dada2)
library(DECIPHER)
library(phangorn)
library(ggplot2)
library(phyloseq)
set.seed(100)

list.files("./")
miseq_path <- "./" 



fnFs <- sort(list.files(miseq_path, pattern=".fastq"))

sampleNames <- sapply(strsplit(fnFs, ".fastq"), `[`, 1)
sampleNames <- sapply(strsplit(sampleNames, "-"), `[`, 2)
head(sampleNames)




sampleNames <- sapply(strsplit(fnFs, "_L001_R1_001.fastq.gz"), `[`,1)
#sampleNames <- sapply(strsplit(basename(fnFs), ".."), `[`, 0)
head(sampleNames)

filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_filt.fastq.gz"))

#filtFs <- sort(list.files(filt_path, pattern=".fastq.gz"))
#filtFs <- file.path(filt_path, paste0(fnFs, "_F_filt.fastq.gz"))
#filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))
# out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
#                      maxN = 0, maxEE = c(2, 2), 
#                      truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread=FALSE,verbose=TRUE) # On Windows set multithread=FALSE
# head(out)


out <- filterAndTrim(fnFs, filtFs,
                     maxN = 0, maxEE = 2, 
                     truncQ = 2, truncLen=80,minLen = 50, rm.phix = TRUE, compress = TRUE, multithread=FALSE,verbose=TRUE) # On Windows set multithread=FALSE


out <- filterAndTrim(fnFs, filtFs,
                     maxN = 0, maxEE = 2, 
                     truncQ = 2, truncLen=100,minLen = 50, compress = TRUE, multithread=FALSE,verbose=TRUE)
head(out)


out <- filterAndTrim(fnFs, filtFs,
                     maxN = 0, maxEE = 2, 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread=FALSE,verbose=TRUE) # On Windows set multithread=FALSE


filterAndTrim(testFastqs, filtFastqs, maxN=0, maxEE=2, verbose=TRUE)
filterAndTrim(testFastqs, filtFastqs, truncQ=2, truncLen=200, rm.phix=TRUE, rm.lowcomplex=8)

filtered=out
plotQualityProfile(fnFs[1])
fnFs[1]
ExtPosU1_S27_L001_R2_001_filt.fastq.gz
filtFs[filtFs != "ExtPosU1_S27_L001_R2_001_filt.fastq.gz"]


errF <- learnErrors(filtFs, multithread=FALSE)
dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
#seqtab <- makeSequenceTable(mergers)
seqtab <- makeSequenceTable(dadaFs)
saveRDS(seqtab, file = "seqtab.rds")

taxa1 <- assignTaxonomy(seqtab, 'D:/base/phil6/files6/silva_nr99_v138.1_train_set.fa.gz', multithread = FALSE, tryRC = TRUE)
taxa2 <- assignTaxonomy(seqtab, 'silva_nr99_v138.1_wSpecies_train_set.fa.gz', multithread = FALSE, tryRC = TRUE)
dev.new()
plotQualityProfile(fnFs[1:2])
getwd()
write.csv(otu_table(seqtab, taxa_are_rows=FALSE),'taxa.csv')

tt.plus <- addSpecies(taxa1, "D:/base/phil6/files6/silva_species_assignment_v138.1.fa.gz", verbose=TRUE)
write.csv(tax_table(tt.plus),'lineage1.csv')



taxa <- assignTaxonomy(seqtab, 'D:/base/sarraj/test/sh_general_release_dynamic_all_29.11.2022.fasta', multithread = TRUE, tryRC = TRUE)
write.csv(otu_table(seqtab, taxa_are_rows=FALSE),'taxa.csv')
write.csv(tax_table(taxa1),'lineage.csv')
