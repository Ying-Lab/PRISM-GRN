## preprocessin
## Convert scATAC-seq data into gene activity score

library(Signac)
library(Seurat)
library(Matrix)
library(data.table)
library(BiocGenerics)
library(GenomicRanges)
library(future) ## 是否调用多线程

library(EnsDb.Hsapiens.v75)
library(BSgenome.Hsapiens.UCSC.hg19)
seq.levels = c(1:22, "X", "Y")
include.body = TRUE
upstream = 2000
downstream = 0
verbose = TRUE

# read scATAC-seq cell*peak sparse matrix
peak_data <- fread("./GSM4156592_GM12878.rep3.counts.txt",skip = 1)
peak.matrix <- sparseMatrix(i = peak_data$V1, 
                              j = peak_data$V2, 
                              x = peak_data$V3)
# read metadata
metacell = read.table("./GSM4156592_GM12878.rep3.barcodes.txt")
metapeak = read.table("./GSM4156592_GM12878.rep3.peaks.bed")

## assign row and col names
colnames(peak.matrix) = metacell$V1
peak_list <- list()
peak_list <- apply(metapeak, 1, function(row) {
  paste(row[1], ":", row[2], "-", row[3])
})
rownames(peak.matrix) = peak_list
colnames(metapeak) = c('chr','start','end')
peak.gr = makeGRangesFromDataFrame(metapeak)

# reference genome assembly
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
start(peak.gr[BiocGenerics::start(peak.gr) == 0, ]) <- 1
# choose the chromatin (depends on species)
seq.levels = c(1:22, "X", "Y") 
gtf <- GenomeInfoDb::keepSeqlevels(x = annotation, value = seq.levels, pruning.mode = 'coarse')

if (!any(GenomeInfoDb::seqlevelsStyle(x = gtf) == GenomeInfoDb::seqlevelsStyle(x = peak.gr))) {
    GenomeInfoDb::seqlevelsStyle(gtf) <- GenomeInfoDb::seqlevelsStyle(peak.gr)
}
gtf.genes <- gtf 

# the genomic distance 
upstream <- 2000
downstream <- 2000
include.body <- TRUE
if (include.body) {
    gtf.body_prom <- Signac::Extend(x = gtf.genes, upstream = upstream, downstream = downstream)
} else {
    gtf.body_prom <- SummarizedExperiment::promoters(x = gtf.genes, upstream = upstream, downstream = downstream)
}

gene.distances <- GenomicRanges::distanceToNearest(x = peak.gr, subject = gtf.body_prom)
keep.overlaps <- gene.distances[rtracklayer::mcols(x = gene.distances)$distance == 0]  
peak.ids <- peak.gr[S4Vectors::queryHits(x = keep.overlaps)]  
gene.ids <- gtf.genes[S4Vectors::subjectHits(x = keep.overlaps)]  
gene.ids$gene_name[is.na(gene.ids$gene_name)] <- gene.ids$gene_id[is.na(gene.ids$gene_name)]
peak.ids$gene.name <- gene.ids$gene_name
peak.ids <- as.data.frame(x = peak.ids)
peak.ids$peak <- rownames(peak.matrix)[S4Vectors::queryHits(x = keep.overlaps)]
annotations <- peak.ids[, c('peak', 'gene.name')]
colnames(x = annotations) <- c('feature', 'new_feature')  

peak.matrix <- as(object = peak.matrix, Class = 'matrix')
all.features <- unique(x = annotations$new_feature)


plan(multisession, workers = 1) 
if (future::nbrOfWorkers() > 1) {
    mysapply <- future.apply::future_sapply
} else {
    mysapply <- ifelse(test = verbose, yes = pbapply::pbsapply, no = sapply)
}


newmat <- mysapply(X = 1:length(x = all.features), FUN = function(x) {
        features.use <- annotations[annotations$new_feature == all.features[[x]], ]$feature
        submat <- peak.matrix[features.use, ]

        if (length(x = features.use) > 1) {
            return(Matrix::colSums(x = submat))  
        } else {
            return(submat)  
        }
    })

# save converted cell*gene TF activity score matrix
newmat <- t(x = newmat)
rownames(x = newmat) <- all.features
colnames(x = newmat) <- colnames(x = peak.matrix)
genescore <- as(object = newmat, Class = 'dgCMatrix')
write.csv(genescore, file = "./gene_activity_rep3.csv", row.names = TRUE,col.names = TRUE)
