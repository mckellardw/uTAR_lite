
args = (commandArgs(TRUE))
if (length(args) == 0) {
  stop("Please specify gene and TAR expression matrices.")
} else {
  TARFile <- args[1] # second argument is TAR File
} 

if (!require('Seurat', quietly=T)) {
  install.packages('Seurat', repo="https://cloud.r-project.org")
}
suppressPackageStartupMessages(library(Seurat,  verbose=F)) # Please use Seurat >= v4.0

if (!require('data.table', quietly=T)) {
  install.packages("data.table", repo="https://cloud.r-project.org")
}
suppressPackageStartupMessages(library(data.table,  verbose=F))


# Borrowed from the Marioni Lab, DropletUtils package (https://rdrr.io/github/MarioniLab/DropletUtils/src/R/write10xCounts.R)
#   (Had R version issues getting it to work as a dependency)
#' @importFrom utils write.table
#' @importFrom Matrix writeMM
#' @importFrom R.utils gzip
.write_sparse <- function(path, x, barcodes, gene.id, gene.symbol, gene.type) {
  # dir.create(path, showWarnings=FALSE)
  gene.info <- data.frame(gene.id, gene.symbol, stringsAsFactors=FALSE)

  gene.info$gene.type <- rep(gene.type, length.out=nrow(gene.info))
  mhandle <- file.path(path, "matrix.mtx")
  bhandle <- gzfile(file.path(path, "barcodes.tsv.gz"), open="wb")
  fhandle <- gzfile(file.path(path, "features.tsv.gz"), open="wb")
  on.exit({
    close(bhandle)
    close(fhandle)
  })

  writeMM(x, file=mhandle)
  write(barcodes, file=bhandle)
  write.table(gene.info, file=fhandle, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

  # Annoyingly, writeMM doesn't take connection objects.
  gzip(mhandle)

  return(NULL)
}



# Read in uTAR matrix, and write to .mtx file
# load TAR matrix and combine with gene expression matrix and subset for valid cells ####
cat("Loading in TAR matrix... ")
TARMat <- fread(
      TARFile,
      sep = "\t",
      header = T,
      stringsAsFactors = F,
      showProgress = F
)

rownames(TARMat) <- tolower(rownames(TARMat))
TARMat <- as.data.frame(TARMat) # force to data frame
rownames(TARMat) <- TARMat$GENE # set the rownames as GENEs
TARMat <- TARMat[, -1] # take out first column
TARMat <- as.sparse(TARMat) #force to sparse to match with gene matrix
cat("Done.\n")
cat("Loaded in ", ncol(TARMat), "cells and ", nrow(TARMat), "TARs.\n")

cat("Saving TAR matrix in MTX format...")
newMatDirName <- paste0(
    stringr::str_remove(TARFile,"TAR_expression_matrix_withDir.tsv.gz"),
    "TAR_feature_bc_matrix"
)
if(!dir.exists(newMatDirName)){
    dir.create(newMatDirName)
}
if(!file.exists(paste0(newMatDirName,"/matrix.mtx.gz"))){
    .write_sparse(
    path=newMatDirName,
    x=TARMat,
    barcodes=colnames(TARMat),
    gene.id=rownames(TARMat),
    gene.symbol=rownames(TARMat),
    gene.type="TAR"
    )
}
cat("Done.\n")