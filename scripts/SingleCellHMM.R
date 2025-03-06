# R --vanilla --slave --args --input_file=PREFIX.bed --output_dir=$(pwd) < SingleCellHMM.R

# Arguments here
args = commandArgs(TRUE)
input_f <- NULL
output_dir <- NULL

for (arg in args) {
  if (grepl("--input_file=", arg)) {
    input_f <- sub("--input_file=", "", arg)
  } else if (grepl("--output_dir=", arg)) {
    output_dir <- sub("--output_dir=", "", arg)
  }
}

if (is.null(input_f) || is.null(output_dir)) {
  stop("Both --input_file and --output_dir arguments must be provided.")
}

setwd(output_dir)

suppressPackageStartupMessages(library(groHMM, quietly = TRUE))
suppressPackageStartupMessages(library(rtracklayer, quietly = TRUE))

# Function to convert GRanges to BED format
GRangeTobed <- function(gr, f_name) {
  df <- data.frame(
    seqnames = seqnames(gr),
    starts = start(gr) - 1,
    ends = end(gr),
    names = c(rep(".", length(gr))),
    scores = c(rep(".", length(gr))),
    strands = strand(gr)
  )
  write.table(
    df,
    file = f_name,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
}

# Function to detect transcripts using HMM
detectTranscripts_AllEM <- function(reads = NULL,
                                    Fp = NULL,
                                    Fm = NULL,
                                    LtProbA = -5,
                                    LtProbB = -200,
                                    UTS = 5,
                                    size = 50,
                                    threshold = 0.1,
                                    debug = TRUE,
                                    ...) {
  stopifnot(!is.null(reads) | (!is.null(Fp) & !is.null(Fm)))
  
  # Setup/Window Analysis/Casting
  epsilon <- 0.001
  
  # Allow equivalent form of Fp and Fm to be specified in the function automatically
  if (is.null(Fp) & is.null(Fm)) {
    Fp <-
      windowAnalysis(reads = reads,
                     strand = "+",
                     windowSize = size,
                     ...)
    Fm <-
      windowAnalysis(reads = reads,
                     strand = "-",
                     windowSize = size,
                     ...)
  }
  
  nFp <- NROW(Fp)
  nFm <- NROW(Fm)
  CHRp <- as.character(names(Fp))
  CHRm <- as.character(names(Fm))
  
  size <- as.integer(size)
  ANS <- NULL
  
  # Set up initial HMM variables
  HMM <- list()
  HMM$nstates <- as.integer(2)
  HMM$ePrDist <- c("dgamma", "dgamma")
  HMM$iProb <-
    as.double(log(c(1.0, 0.0))) # Non-transcribed, transcribed
  HMM$ePrVars <-
    as.list(data.frame(c(UTS, 1 / UTS,-1), c(0.5, 10,-1)))
  HMM$tProb <- as.list(data.frame(c(log(1 - exp(
    LtProbA
  )), LtProbA),
  c(LtProbB, log(1 - exp(
    LtProbB
  )))))
  
  # Cast counts to a real, and combine +/- strand into one list variable
  FT <- list()
  for (i in 1:nFp)
    FT[[i]] <- as.double(Fp[[i]] + 1)
  for (i in 1:nFm)
    FT[[i + nFp]] <- as.double(Fm[[i]] + 1)
  
  # Free unused memory
  remove(Fp)
  remove(Fm)
  
  # Run EM algorithm
  BWem <- .Call(
    "RBaumWelchEM",
    HMM$nstates,
    FT,
    as.integer(1),
    HMM$ePrDist,
    HMM$ePrVars,
    HMM$tProb,
    HMM$iProb,
    as.double(threshold),
    c(TRUE, TRUE),
    c(TRUE, TRUE),
    as.integer(1),
    TRUE,
    PACKAGE = "groHMM"
  )
  
  # Translate these into transcript positions
  for (i in seq_along(CHRp)) {
    ans <- .Call("getTranscriptPositions",
                 as.double(BWem[[3]][[i]]),
                 as.double(0.5),
                 size,
                 PACKAGE = "groHMM")
    Nrep <- NROW(ans$Start)
    ANS <- rbind(ANS,
                 data.frame(
                   chrom = rep(CHRp[i], Nrep),
                   start = ans$Start,
                   end = ans$End,
                   strand = rep("+", Nrep)
                 ))
  }
  
  for (i in seq_along(CHRm)) {
    ans <-
      .Call("getTranscriptPositions",
            as.double(BWem[[3]][[i + nFp]]),
            as.double(0.5),
            size,
            PACKAGE = "groHMM")
    Nrep <- NROW(ans$Start)
    ANS <-
      rbind(ANS,
            data.frame(
              chrom = rep(CHRm[i], NROW(ans$Start)),
              start = ans$Start,
              end = ans$End,
              strand = rep("-", Nrep)
            ))
  }
  
  BWem[[4]] <- GRanges(
    seqnames = Rle(ANS$chrom),
    ranges = IRanges(ANS$start, ANS$end - 1),
    strand = Rle(strand(ANS$strand)),
    type = Rle("tx", NROW(ANS)),
    ID = paste(ANS$chrom, "_", ANS$start, ANS$strand, sep = "")
  )
  names(BWem) <-
    c("emisParams", "transParams", "viterbiStates", "transcripts")
  
  if (debug) {
    print(BWem[[1]])
    print(BWem[[2]])
  }
  
  return(BWem)
}

# Main script execution
S_split = import(input_f)
hmmResult_AllEM_split <- detectTranscripts_AllEM(S_split)
output_file <-
  file.path(output_dir, paste(strsplit(basename(input_f), '.bed', TRUE)[[1]][1], "_HMM.bed", sep = ""))
GRangeTobed(hmmResult_AllEM_split$transcripts, output_file)

cat("Checking for output .bed file from R script... ")
cat(file.exists(output_file), "\n")
