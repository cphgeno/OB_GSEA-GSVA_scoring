suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(pheatmap)
  library(GSVA)
  library(qusage)
})

# source helper functions
get_script_dir <- function() {
  commandArgs() %>%
    tibble::enframe(name = NULL) %>%
    tidyr::separate(
      col = value, into = c("key", "value"), sep = "=", fill = "right"
    ) %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value) %>%
    dirname(.)
}
source(file.path(get_script_dir(), "gsva_helper_functions.R"))

# parse options
option_list <- list(
  make_option("--output_dir", "-o",
    dest="output_dir", type="character",
    help="output directory where files will be saved"),
  make_option("--name",
    type = "character",
    help = "name of the dataset"),
  make_option("--preprocessing.sampleranks",
    dest="ranks_df", type="character",
    help="input directory where datasets and genesets are found"),
  make_option("--preprocessing.samplediffranks",
    dest="ranksdiff_df", type="character",
    help="input directory where datasets and genesets are found"),
  make_option("--preprocessing.referenceranks",
    dest="reference_df", type="character",
    help="input directory where datasets and genesets are found"),
  make_option("--preprocessing.meta",
    dest="metadata", type="character",
    help="input directory where datasets and genesets are found"),
  make_option("--preprocessing.genesets",
    dest="genesets", type="character",
    help="input directory where datasets and genesets are found"),
  make_option("--input_type",
    dest="input_type", type="character",
    help="input_type depending on ranked input"),
  make_option("--algorithm",
    type = "character", default = NULL,
    help = "Algorithm to be used for analysis")
)
opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)


# --- Required argument check ---
defined_opts <- sapply(option_list, function(x) x@dest)
missing_opts <- defined_opts[vapply(opts[defined_opts], is.null, logical(1))]
if (length(missing_opts) > 0) {
  print_help(opt_parser)
  stop(paste("Missing required option(s):",
             paste(missing_opts, collapse=", ")),
       call.=FALSE)
}

TOOL_COLOURS <- list(
  gsva_RankReference = "#05A005",
  plage_RankReference = "#E93CB5",
  zscore_RankReference = "#FF7E0D",
  gsva_RankExpr = "#47B847",
  gsva_DeltaCentroid  = "#237F23",
  plage_RankExpr = "#E691D0",
  plage_DeltaCentroid  = "#B55F9B",
  zscore_RankExpr = "#FF993E",
  zscore_DeltaCentroid  = "#CC6624"
)

tool_colour = TOOL_COLOURS[[paste0(opts$algorithm, "_", opts$input_type)]]

metadata_df <- read.table(opts$metadata, sep = '\t', header = TRUE, check.names = FALSE)
genesets_list <- read.gmt(opts$genesets)

# check analysis type to fetch correct input and run correct class wrapper
if (opts$input_type == "RankExpr"){
  input_df <- opts$ranks_df_diff
  gsva_wrapper_classI_III(input_df, metadata_df, genesets_list, opts$name, opts$output_dir, opts$input_type, tool_colour)
} else if (opts$input_type == "DeltaCentroid"){
  input_df <- opts$ranks_df_diff
  gsva_wrapper_classI_III(input_df, metadata_df, genesets_list, opts$name, opts$output_dir, opts$input_type, tool_colour)
} else if (opts$input_type == "RankReference"){
  input_df <- read.table(opts$ranks_df, sep = '\t', header = TRUE, check.names = FALSE)
  input_df <- column_to_rownames(input_df, var = 'Geneid')
  reference_df <- read.table(opts$reference_df, sep = '\t', header = TRUE, check.names = FALSE)
  reference_df <- column_to_rownames(reference_df, var = 'Geneid')
  gsva_wrapper_classII(opts$output_dir, input_df, reference_df, metadata_df, genesets_list, opts$algorithm, opts$name, tool_colour)
} else {
  stop(paste("Input type", opts$input_type, "not allowed for this tool"))
}

