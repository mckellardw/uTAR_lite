#TODO- Add description for this function
checkIfExistGene_noDir <- function(input,gene_ref){
  chrom <- input[[1]]
  startPos <- as.numeric(input[[2]])
  endPos <- as.numeric(input[[3]])
  gene_ref <- gene_ref # bed of entire genes

  # 1: either partially in or completely inside gene
  chrMatch_gene <- (chrom==gene_ref$chr)
  qinrstartMatch <- (startPos>=gene_ref$start)&(startPos<=gene_ref$end)
  qinrendMatch <- (endPos<=gene_ref$end)&(endPos>=gene_ref$start)
  qinroutGeneAll_gene <- (chrMatch_gene*qinrstartMatch*qinrendMatch) # query completely inside ref
  rinqstartMatch <- (startPos<(gene_ref$start))
  rinqendMatch <- (endPos>(gene_ref$end))
  rinqoutGeneAll <- (chrMatch_gene*rinqstartMatch*rinqendMatch) # ref completely inside query

  # also check partial overlap from beginning
  partialStart <- (chrMatch_gene*rinqstartMatch*qinrendMatch)
  partialEnd <- (chrMatch_gene*qinrstartMatch*rinqendMatch)
  if(sum(qinroutGeneAll_gene)|sum(rinqoutGeneAll)|sum(partialStart)|sum(partialEnd)){
    # g's are for gene names
    g1 <- as.character(gene_ref$gene[as.logical(qinroutGeneAll_gene)])
    g2 <- as.character(gene_ref$gene[as.logical(rinqoutGeneAll)])
    g3 <- as.character(gene_ref$gene[as.logical(partialStart)])
    g4 <- as.character(gene_ref$gene[as.logical(partialEnd)])
    gAll <- unique(c(g1,g2,g3,g4))
    return(paste(c(gAll,"1"),collapse="_"))
  }else{
    return(0)
  }
}

#TODO- Add description for this function, and reqrite if possible...
checkIfExistGene2 <- function(
  input,
  gene_ref
){
  chrom <- input[[1]]
  startPos <- as.numeric(input[[2]])
  endPos <- as.numeric(input[[3]])
  direction <- input[[4]]
  gene_ref <- gene_ref # bed of entire genes
  #print(paste(chrom,startPos,endPos))

  # 1: either partially in or completely inside gene
  chrMatch_gene <- (chrom==gene_ref$chr)
  dirMatch_gene <- (direction==gene_ref$direction)
  qinrstartMatch <- (startPos>=gene_ref$start)&(startPos<=gene_ref$end)
  qinrendMatch <- (endPos<=gene_ref$end)&(endPos>=gene_ref$start)
  qinroutGeneAll_gene <- (chrMatch_gene*qinrstartMatch*qinrendMatch*dirMatch_gene) # query completely inside ref
  rinqstartMatch <- (startPos<(gene_ref$start))
  rinqendMatch <- (endPos>(gene_ref$end))
  rinqoutGeneAll <- (chrMatch_gene*rinqstartMatch*rinqendMatch*dirMatch_gene) # ref completely inside query
  # also check partial overlap from beginning
  partialStart <- (chrMatch_gene*rinqstartMatch*qinrendMatch*dirMatch_gene)
  partialEnd <- (chrMatch_gene*qinrstartMatch*rinqendMatch*dirMatch_gene)
  if(sum(qinroutGeneAll_gene)|sum(rinqoutGeneAll)|sum(partialStart)|sum(partialEnd)){
    # g's are for gene names
    g1 <- as.character(gene_ref$gene[as.logical(qinroutGeneAll_gene)])
    g2 <- as.character(gene_ref$gene[as.logical(rinqoutGeneAll)])
    g3 <- as.character(gene_ref$gene[as.logical(partialStart)])
    g4 <- as.character(gene_ref$gene[as.logical(partialEnd)])
    gAll <- unique(c(g1,g2,g3,g4))

    # d's are for directions
    d1 <- as.character(gene_ref$direction[as.logical(qinroutGeneAll_gene)])
    d2 <- as.character(gene_ref$direction[as.logical(rinqoutGeneAll)])
    d3 <- as.character(gene_ref$direction[as.logical(partialStart)])
    d4 <- as.character(gene_ref$direction[as.logical(partialEnd)])
    dAll <- unique(c(d1,d2,d3,d4))

    return(paste(c(gAll,dAll,"1"),collapse="_"))
  }else{
    #return(paste(c(direction,"0"),collapse="_"))
    return(0)
  }
}


# MAKE SURE CHROM IN REFERENCE AND HMM BED MATCH EACH OTHER
#TODO- refactor so we don't need refFlat format, and can directly load in .gtf
args=(commandArgs(TRUE))
inputRefFlat <- args[1]
HMMbedFile <- args[2]
num_cores <- as.numeric(args[3])

refName <- unlist(strsplit(inputRefFlat,"/"))[length(unlist(strsplit(inputRefFlat,"/")))]
# read in ref genes
gene_ref <- read.delim(
  inputRefFlat,
  header=F,
  comment.char="#"
)
gene_ref <- gene_ref[,c("V1","V3","V4","V5","V6")] ####### edit this based on file type
colnames(gene_ref) <- c("gene","chr","direction","start","end")
gene_ref$gene <- as.character(gene_ref$gene)
gene_ref$chr <- as.character(gene_ref$chr)
gene_ref$direction <- as.character(gene_ref$direction)
gene_ref$start <- as.numeric(gene_ref$start)
gene_ref$end <- as.numeric(gene_ref$end)

input <- HMMbedFile
HMManno <- read.delim(input, header=FALSE)
HMManno_bare <- HMManno[,c("V1","V2","V3","V6")]
colnames(HMManno_bare) <- c("chr","start","end","direction")
HMManno_bare$chr <- as.character(HMManno_bare$chr)
HMManno_bare$direction <- as.character(HMManno_bare$direction)
HMManno_bare$start <- as.numeric(HMManno_bare$start)
HMManno_bare$end <- as.numeric(HMManno_bare$end)

########################################################

if(!require("parallel")){
	install.packages("parallel")
}
library(parallel, quietly=T)

clust <- makeCluster(num_cores)

HMManno$inGene <- parApply(
  cl = clust,
  X=HMManno_bare,
  MARGIN=1,
  FUN=checkIfExistGene_noDir,
  gene_ref=gene_ref
)
stopCluster(clust)

#HMManno$inGene <- apply(X=HMManno_bare,MARGIN=1,FUN=checkIfExistGene_noDir,gene_ref=gene_ref)

# make dataframe into refFlat file format
HMManno$geneName <- paste(HMManno$V1,"_",HMManno$V2,"_",HMManno$V3,"_",HMManno$V6,"_",HMManno$V7,"_",HMManno$inGene,sep="")
#HMManno$name <- HMManno$geneName
HMManno$name <- paste(HMManno$V1,"_",HMManno$V2,"_",HMManno$V3,sep="")
HMManno$chrom <- HMManno$V1
HMManno$strand <- HMManno$V6
HMManno$txStart <- HMManno$V2
HMManno$txEnd <- HMManno$V3
HMManno$cdsStart <- HMManno$V2
HMManno$cdsEnd <- HMManno$V3
HMManno$exonCount <- 1
HMManno$exonStarts <- HMManno$V2
HMManno$exonEnds <- HMManno$V3

HMMannoReady <- HMManno[,c("geneName","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds")]

#TODO: export as GTF
outFile=paste0(stringr::str_remove(string = input, pattern=".bed.gz"),".noDir.","refFlat")
cat("Writing unstranded TAR annotations (annotated & unannotated) to ", outFile, "\n")
write.table(
  HMMannoReady,
  outFile,
  sep="\t",
  row.names=F,
  col.names = F,
  quote=F
)

###################################################################################
#### include direction (strandedness) below
###################################################################################

#HMManno_bare_sample <- df[sample(nrow(df),10),]
#HMManno$inAnno <- apply(X=HMManno_bare,MARGIN=1,FUN=checkIfExist,exon_ref=gene_gtf_bare,gene_ref=gene_ref)
#HMManno$inGene <- apply(X=HMManno_bare,MARGIN=1,FUN=checkIfExistGene2,gene_ref=gene_ref)

library(parallel,quietly=T)

clust <- makeCluster(num_cores)
HMManno$inGene <- parApply(
  cl = clust,
  X=HMManno_bare,
  MARGIN=1,
  FUN=checkIfExistGene2,
  gene_ref=gene_ref
)
stopCluster(clust)

# make dataframe into refFlat file format
HMManno$geneName <- paste(HMManno$V1,"_",HMManno$V2,"_",HMManno$V3,"_",HMManno$V6,"_",HMManno$V7,"_",HMManno$inGene,sep="")
#HMManno$name <- HMManno$geneName
HMManno$name <- paste(HMManno$V1,"_",HMManno$V2,"_",HMManno$V3,sep="")
HMManno$chrom <- HMManno$V1
HMManno$strand <- HMManno$V6
HMManno$txStart <- HMManno$V2
HMManno$txEnd <- HMManno$V3
HMManno$cdsStart <- HMManno$V2
HMManno$cdsEnd <- HMManno$V3
HMManno$exonCount <- 1
HMManno$exonStarts <- HMManno$V2
HMManno$exonEnds <- HMManno$V3

HMManno <- HMManno[,c("geneName","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds")]

# Filter to only include unannotated TARs
HMManno <- HMManno[ stringr::str_ends(HMManno$geneName,pattern = "_0"),] # Note- suffix changes to "-0" when loaded into Seurat

if(nrow(HMManno) == 0){
  message("No uTARs found!")
}

#TODO export as GTF
outFile=paste0(stringr::str_remove(string = input, pattern=".bed.gz"),".withDir.","refFlat")
cat("Writing stranded uTAR annotations to ", outFile, "\n")
write.table(
  HMManno,
  outFile,
  sep="\t",
  row.names=F,
  col.names = F,
  quote=F
)
