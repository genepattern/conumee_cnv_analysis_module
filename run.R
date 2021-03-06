write("Loading dependencies...", stdout())
suppressMessages(suppressWarnings(library("optparse")))
suppressMessages(suppressWarnings(library("conumee")))
suppressMessages(suppressWarnings(library("rtracklayer")))
suppressMessages(suppressWarnings(library("foreach")))
suppressMessages(suppressWarnings(library("doParallel")))
suppressMessages(suppressWarnings(library("parallel")))
suppressMessages(suppressWarnings(library("biomaRt")))

# Parse input arguments
parser <- OptionParser()
parser <- add_option(parser, c("--data"), help = "methylation data")
parser <- add_option(parser, c("--controls"), help = "names of control samples")
parser <- add_option(parser, c("--controlsdata"), help = "controls methylation data")
parser <- add_option(parser, c("--arraytype"), help = "450k, EPIC, or both")
parser <- add_option(parser, c("--genesfile"), help = "file with list of genes to highlight")
parser <- add_option(parser, c("--ignorefile"), help = "bed file of regions to exclude")
parser <- add_option(parser, c("--xy"), help = "include XY chromosomes")
args <- parse_args(parser)

# methylation data
write("Loading data...", stdout())

all.data <- readRDS(args$data)
if (class(all.data) != "MethylSet") {
  write("Data class is not a MethylSet.", stdout())
  stop()
}

# Array type parameters
if (args$arraytype == "450k") {
  arraytype <- "450k"
} else if (args$arraytype == "EPIC") {
  arraytype <- "EPIC"
} else if (args$arraytype == "both") {
  arraytype <- "overlap"
} else {
  write("Array type not recognized.", stdout())
  stop()
}

# Control sample names
if (tolower(args$controls) != "none") {
  write("Using subset as controls...", stdout())
  controls.names <- strsplit(args$controls, ",")
  controls.names <- lapply(controls.names, trimws)
  controls.names <- unlist(controls.names)

  # bugfix conumee::CNV.fit cor() step
  if (length(controls.names) == 1) {
    controls.names <- c(controls.names, controls.names)
  }

  if (!all(controls.names %in% colnames(all.data))) {
    write("One or more of the provided names of control samples are not in the dataset.",
      stdout())
    stop()
  }
} else {
  if (is.null(args$controlsdata)) {
    write("Please specify control sample names or provide a separate controls methylation data file.",
      stdout())
    stop()
  }
  # Load control data
  write("Loading separate controls data...", stdout())
  controls.data <- readRDS(args$controlsdata)
  controls.names <- colnames(controls.data)

  if (arraytype == "450k" || arraytype == "overlap") {
    all.data <- combineArrays(all.data, controls.data, outType = "IlluminaHumanMethylation450k")
  } else {
    all.data <- combineArrays(all.data, controls.data, outType = "IlluminaHumanMethylationEPIC")
  }
}

# Query biomart for gene annotations
if (!("genesfile" %in% names(args))) {
  detail_regions <- NULL
} else {
  write("Querying biomaRt for gene locations...", stdout())
  gene.list <- read.table(args$genesfile, header = FALSE)
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "grch37.ensembl.org")
  detail_regions <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position",
    "end_position", "strand"), filter = "hgnc_symbol", values = gene.list, mart = ensembl)

  valid_regions <- detail_regions$chromosome_name %in% seq(1, 22)
  detail_regions <- detail_regions[valid_regions, ]

  # Include 5kb upstream for promoter region
  for (i in 1:dim(detail_regions)[1]) {
    entry <- detail_regions[i, ]
    if (entry$strand == 1) {
      entry$start_position <- entry$start_position - 5000
    } else {
      entry$end_position <- entry$end_position + 5000
    }
  }

  # Convert to GRange object
  detail_regions$strand[detail_regions$strand == 1] <- "+"
  detail_regions$strand[detail_regions$strand == -1] <- "-"
  detail_regions$chromosome_name <- paste("chr", detail_regions$chromosome_name,
    sep = "")
  detail_regions <- makeGRangesFromDataFrame(detail_regions, seqnames = "chromosome_name",
    start.field = "start_position", end.field = "end_position", strand.field = "strand",
    keep.extra.columns = TRUE)
  detail_regions$name <- detail_regions$hgnc_symbol
  detail_regions$hgnc_symbol <- NULL
}

# Include/exclude sex chromosomes
if (tolower(args$xy) == "yes") {
  xy <- TRUE
} else {
  xy <- FALSE
}

# Ignore highly polymorphic regions of genome
if (!("ignorefile" %in% names(args))) {
  ignore_regions <- NULL
} else if (endsWith(args$ignorefile, ".bed")) {
  write("Loading regions to exclude...", stdout())
  ignore_regions <- args$ignorefile
} else {
  write(paste("Ignore regions value, '", args$ignorefile, "' not recognized.",
    sep = ""), stdout())
  stop()
}

# Create annotation object. Set to look at probes common to both 850k and 450k.
write("Creating CNV annotation object...", stdout())
anno <- CNV.create_anno(array_type = arraytype, exclude_regions = ignore_regions,
  detail_regions = detail_regions, chrXY = xy, bin_minprobes = 15, bin_minsize = 50000,
  bin_maxsize = 5e+06)

# Create CNV object from methylation data
write("Creating CNV object...", stdout())

# Create CNV object
cnv.data <- CNV.load(all.data)

# parallelization parameters TODO use max cores
all_samples <- colnames(cnv.data@intensity)
numCores <- min(length(all_samples) - length(controls.names), detectCores() - 2)
cl <- makeCluster(numCores)
registerDoParallel(cl)

# CNV analysis/plot function
cnv.analyze.plot <- function(sample, controls.names, cnv.data, anno) {
  write(paste(sample, " - fitting", sep = ""), stdout())
  cnv.analysis <- CNV.fit(cnv.data[sample], cnv.data[controls.names], anno)
  cnv.analysis <- CNV.bin(cnv.analysis)
  cnv.analysis <- CNV.detail(cnv.analysis)
  cnv.analysis <- CNV.segment(cnv.analysis)
  CNV.detail(cnv.analysis)
  # Open file to write
  write(paste(sample, " - plotting", sep = ""), stdout())
  plot.filename <- paste(sample, ".cnvPlots.pdf", sep = "")

  pdf(plot.filename, width = 16, height = 8)

  # Overview plots
  CNV.genomeplot(cnv.analysis)
  CNV.detailplot_wrap(cnv.analysis)

  # Chromosome plots
  for (chrom in rownames(anno@genome)) {
    CNV.genomeplot(cnv.analysis, chrom)
  }

  # Detail plots
  genes <- names(cnv.analysis@detail$probes)
  for (gene in genes) {
    CNV.detailplot(cnv.analysis, gene)
  }
  dev.off()

  # Export seg file
  CNV.write(cnv.analysis, what = "segments", file = paste(sample, ".cnv.seg", sep = ""))

  # Export details file
  CNV.write(cnv.analysis, what = "detail", file = paste(sample, ".detail.cnv.seg",
    sep = ""))
}

# Analyze in parallel
write(paste("Using ", numCores, " cores...", sep = ""), stdout())
foreach(sample = all_samples, .packages = c("conumee")) %dopar% {
  if (!(sample %in% controls.names)) {
    cnv.analyze.plot(sample, controls.names, cnv.data, anno)
  }
}

write("Done.", stdout())
