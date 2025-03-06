

library(parallel, quietly = TRUE)
library(readr, quietly = TRUE)


# Function to check if a gene exists without considering direction
# Args:
#   input: A list containing chromosome, start position, and end position
#   gene_ref: A data frame containing gene reference information
# Returns:
#   A string with gene names if a match is found, otherwise 0
checkIfExistGene_noDir <- function(input, gene_ref) {
  chrom <- input[[1]]
  startPos <- as.numeric(input[[2]])
  endPos <- as.numeric(input[[3]])
  gene_ref <- gene_ref # bed of entire genes
  
  # 1: either partially in or completely inside gene
  chrMatch_gene <- (chrom == gene_ref$chr)
  qinrstartMatch <-
    (startPos >= gene_ref$start) & (startPos <= gene_ref$end)
  qinrendMatch <-
    (endPos <= gene_ref$end) & (endPos >= gene_ref$start)
  qinroutGeneAll_gene <-
    (chrMatch_gene * qinrstartMatch * qinrendMatch) # query completely inside ref
  rinqstartMatch <- (startPos < gene_ref$start)
  rinqendMatch <- (endPos > gene_ref$end)
  rinqoutGeneAll <-
    (chrMatch_gene * rinqstartMatch * rinqendMatch) # ref completely inside query
  
  # also check partial overlap from beginning
  partialStart <- (chrMatch_gene * rinqstartMatch * qinrendMatch)
  partialEnd <- (chrMatch_gene * qinrstartMatch * rinqendMatch)
  if (sum(qinroutGeneAll_gene) |
      sum(rinqoutGeneAll) | sum(partialStart) | sum(partialEnd)) {
    # g's are for gene names
    g1 <-
      as.character(gene_ref$gene[as.logical(qinroutGeneAll_gene)])
    g2 <- as.character(gene_ref$gene[as.logical(rinqoutGeneAll)])
    g3 <- as.character(gene_ref$gene[as.logical(partialStart)])
    g4 <- as.character(gene_ref$gene[as.logical(partialEnd)])
    gAll <- unique(c(g1, g2, g3, g4))
    return(paste(c(gAll, "1"), collapse = "_"))
  } else {
    return(0)
  }
}

# Function to check if a gene exists considering direction
# Args:
#   input: A list containing chromosome, start position, end position, and direction
#   gene_ref: A data frame containing gene reference information
# Returns:
#   A string with gene names and directions if a match is found, otherwise 0
checkIfExistGene2 <- function(input, gene_ref) {
  chrom <- input[[1]]
  startPos <- as.numeric(input[[2]])
  endPos <- as.numeric(input[[3]])
  direction <- input[[4]]
  gene_ref <- gene_ref # bed of entire genes
  
  # 1: either partially in or completely inside gene
  chrMatch_gene <- (chrom == gene_ref$chr)
  dirMatch_gene <- (direction == gene_ref$direction)
  qinrstartMatch <-
    (startPos >= gene_ref$start) & (startPos <= gene_ref$end)
  qinrendMatch <-
    (endPos <= gene_ref$end) & (endPos >= gene_ref$start)
  qinroutGeneAll_gene <-
    (chrMatch_gene * qinrstartMatch * qinrendMatch * dirMatch_gene) # query completely inside ref
  rinqstartMatch <- (startPos < gene_ref$start)
  rinqendMatch <- (endPos > gene_ref$end)
  rinqoutGeneAll <-
    (chrMatch_gene * rinqstartMatch * rinqendMatch * dirMatch_gene) # ref completely inside query
  
  # also check partial overlap from beginning
  partialStart <-
    (chrMatch_gene * rinqstartMatch * qinrendMatch * dirMatch_gene)
  partialEnd <-
    (chrMatch_gene * qinrstartMatch * rinqendMatch * dirMatch_gene)
  if (sum(qinroutGeneAll_gene) |
      sum(rinqoutGeneAll) | sum(partialStart) | sum(partialEnd)) {
    # g's are for gene names
    g1 <-
      as.character(gene_ref$gene[as.logical(qinroutGeneAll_gene)])
    g2 <- as.character(gene_ref$gene[as.logical(rinqoutGeneAll)])
    g3 <- as.character(gene_ref$gene[as.logical(partialStart)])
    g4 <- as.character(gene_ref$gene[as.logical(partialEnd)])
    gAll <- unique(c(g1, g2, g3, g4))
    
    # d's are for directions
    d1 <-
      as.character(gene_ref$direction[as.logical(qinroutGeneAll_gene)])
    d2 <-
      as.character(gene_ref$direction[as.logical(rinqoutGeneAll)])
    d3 <- as.character(gene_ref$direction[as.logical(partialStart)])
    d4 <- as.character(gene_ref$direction[as.logical(partialEnd)])
    dAll <- unique(c(d1, d2, d3, d4))
    
    return(paste(c(gAll, dAll, "1"), collapse = "_"))
  } else {
    return(0)
  }
}

# MAKE SURE CHROM IN REFERENCE AND HMM BED MATCH EACH OTHER
args <- commandArgs(TRUE)
inputGTF <- args[1]
HMMbedFile <- args[2]
num_cores <- as.numeric(args[3])

# Function to read GTF file and extract relevant columns
# Args:
#   gtf_file: Path to the GTF file
# Returns:
#   A data frame with columns: gene, chr, direction, start, end
readGTF <- function(gtf_file) {
  gtf <-
    readr::read_delim(gtf_file,
                      delim = "\t",
                      comment = "#",
                      col_names = FALSE)
  gtf <- gtf[gtf$X3 == "gene",] # Filter for gene entries
  gtf <- gtf[, c("X9", "X1", "X7", "X4", "X5")]
  colnames(gtf) <- c("gene", "chr", "direction", "start", "end")
  gtf$gene <-
    sapply(strsplit(as.character(gtf$gene), ";"), function(x)
      sub("gene_name ", "", x[grep("gene_name", x)]))
  gtf$chr <- as.character(gtf$chr)
  gtf$direction <- as.character(gtf$direction)
  gtf$start <- as.numeric(gtf$start)
  gtf$end <- as.numeric(gtf$end)
  return(gtf)
}

# Read in gene reference from GTF file
gene_ref <- readGTF(inputGTF)

input <- HMMbedFile
HMManno <- readr::read_delim(input, delim = "\t", col_names = FALSE)
HMManno_bare <- HMManno[, c("X1", "X2", "X3", "X6")]
colnames(HMManno_bare) <- c("chr", "start", "end", "direction")
HMManno_bare$chr <- as.character(HMManno_bare$chr)
HMManno_bare$direction <- as.character(HMManno_bare$direction)
HMManno_bare$start <- as.numeric(HMManno_bare$start)
HMManno_bare$end <- as.numeric(HMManno_bare$end)


clust <- makeCluster(num_cores)

HMManno$inGene <- parApply(
  cl = clust,
  X = HMManno_bare,
  MARGIN = 1,
  FUN = checkIfExistGene_noDir,
  gene_ref = gene_ref
)
stopCluster(clust)

# make dataframe into refFlat file format
HMManno$geneName <-
  paste(
    HMManno$X1,
    "_",
    HMManno$X2,
    "_",
    HMManno$X3,
    "_",
    HMManno$X6,
    "_",
    HMManno$X7,
    "_",
    HMManno$inGene,
    sep = ""
  )
HMManno$name <-
  paste(HMManno$X1, "_", HMManno$X2, "_", HMManno$X3, sep = "")
HMManno$chrom <- HMManno$X1
HMManno$strand <- HMManno$X6
HMManno$txStart <- HMManno$X2
HMManno$txEnd <- HMManno$X3
HMManno$cdsStart <- HMManno$X2
HMManno$cdsEnd <- HMManno$X3
HMManno$exonCount <- 1
HMManno$exonStarts <- HMManno$X2
HMManno$exonEnds <- HMManno$X3

HMMannoReady <-
  HMManno[, c(
    "geneName",
    "name",
    "chrom",
    "strand",
    "txStart",
    "txEnd",
    "cdsStart",
    "cdsEnd",
    "exonCount",
    "exonStarts",
    "exonEnds"
  )]

#TODO: export as GTF
outFile <-
  paste0(stringr::str_remove(string = input, pattern = ".bed.gz"),
         ".noDir.",
         "refFlat")
cat("Writing unstranded TAR annotations (annotated & unannotated) to ",
    outFile,
    "\n")
readr::write_delim(
  HMMannoReady,
  outFile,
  delim = "\t",
  col_names = FALSE,
  quote = FALSE
)


clust <- makeCluster(num_cores)
HMManno$inGene <- parApply(
  cl = clust,
  X = HMManno_bare,
  MARGIN = 1,
  FUN = checkIfExistGene2,
  gene_ref = gene_ref
)
stopCluster(clust)

# make dataframe into refFlat file format
HMManno$geneName <-
  paste(
    HMManno$X1,
    "_",
    HMManno$X2,
    "_",
    HMManno$X3,
    "_",
    HMManno$X6,
    "_",
    HMManno$X7,
    "_",
    HMManno$inGene,
    sep = ""
  )
HMManno$name <-
  paste(HMManno$X1, "_", HMManno$X2, "_", HMManno$X3, sep = "")
HMManno$chrom <- HMManno$X1
HMManno$strand <- HMManno$X6
HMManno$txStart <- HMManno$X2
HMManno$txEnd <- HMManno$X3
HMManno$cdsStart <- HMManno$X2
HMManno$cdsEnd <- HMManno$X3
HMManno$exonCount <- 1
HMManno$exonStarts <- HMManno$X2
HMManno$exonEnds <- HMManno$X3

HMManno <-
  HMManno[, c(
    "geneName",
    "name",
    "chrom",
    "strand",
    "txStart",
    "txEnd",
    "cdsStart",
    "cdsEnd",
    "exonCount",
    "exonStarts",
    "exonEnds"
  )]

# Filter to only include unannotated TARs
HMManno <-
  HMManno[stringr::str_ends(HMManno$geneName, pattern = "_0"), ] # Note- suffix changes to "-0" when loaded into Seurat

if (nrow(HMManno) == 0) {
  message("No uTARs found!")
}

#TODO export as GTF
outFile <-
  paste0(stringr::str_remove(string = input, pattern = ".bed.gz"),
         ".withDir.",
         "refFlat")
cat("Writing stranded uTAR annotations to ", outFile, "\n")
readr::write_delim(
  HMManno,
  outFile,
  delim = "\t",
  col_names = FALSE,
  quote = FALSE
)
