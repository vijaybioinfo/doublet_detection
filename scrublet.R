#!/usr/bin/R

#####################
# Doublet detection #
#####################

# This code will call a method to detect doublets.
# conda activate doublets

library(optparse)
#library(reticulate)
#scr = reticulate::import("scrublet")

optlist <- list(
  make_option(
    opt_str = c("-y", "--yaml"), type = "character", default = "config.yaml",
    help = "Configuration file."
  ),
  make_option(
    opt_str = c("-l", "--log"), default = TRUE,
    help = "Log file: Create a log file rather to stdout."
  ),
  make_option(
    opt_str = c("-v", "--verbose"), type = "logical", default = TRUE,
    help = "Verbose."
  )
)
optparse <- OptionParser(option_list = optlist)
opt <- parse_args(optparse)

# If interactive
#opt$yaml <- "/home/kmlanderos/mesothelioma/scripts/doubdet_tcells.yaml"

if(requireNamespace("crayon", quietly = TRUE)){
  cyan = crayon::cyan; red_bold = crayon::red$bold
}else{ cyan = red_bold = c }

if(opt$verbose){
  cat(red_bold("-----------------------------------------\n"))
  cat(red_bold("----------- Doublet detection -----------\n"))
  cat(red_bold("-----------------------------------------\n"))
}

config = yaml::read_yaml(opt$yaml)

if(opt$verbose) cat(cyan('\n----------------------- Digesting parameters\n'))
parameters2set = c(
  "filters$score" = "0.2", "filters$doublet_rate" = "0.06",
  "filters$lib" = "'Gex'", "vars$class" = "'ht_classification.global'",
  "do_parallel" = "TRUE"
)
for(i in 1:length(parameters2set)){
  param = paste0("config$", names(parameters2set)[i])
  command = paste0(param, " = ", parameters2set[[i]])
  if(is.null(eval(parse(text = param)))){
    if(opt$verbose) cat(" ", command, "\n"); eval(parse(text = command))
  }
}

if(opt$verbose) cat(cyan('\n----------------------- Sourcing code\n'))
packages_funcs = c(
  "ggplot2", "cowplot", "patchwork", "foreach",
  "/home/ciro/scripts/handy_functions/devel/file_reading.R", # readfile
  #"/home/ciro/scripts/doubdet/scrublet_fun.R", #doublets_scrublet, #doublets_scrublet_plot
  "/home/kmlanderos/scripts/doubdet/scrublet_fun.R", #doublets_scrublet, #doublets_scrublet_plot
  "/home/ciro/scripts/clustering/R/utilities.R" # get_source_data
)
loaded <- lapply(X = packages_funcs, FUN = function(x){
    cat("*", x, "\n")
    if(!file.exists(x)){
      suppressMessages(require(package = x, quietly = TRUE, character.only = TRUE))
    }else{ source(x) }
}); theme_set(theme_cowplot())
get_library_id = function(mdata){
  y = c("origlib", "library")
  tvar <- y %in% colnames(mdata)
  if(sum(tvar) != 0) y[tvar]
}

if(opt$verbose) cat(cyan('\n----------------------- Directories structure\n'))
if(basename(config$output_dir) != config$project_name)
  config$output_dir = paste0(config$output_dir, "/", config$project_name)
if(!grepl("scratch|beegfs", getwd())){
  cat("No scratch folder involved; careful about temp files...\n")
  dir.create(config$output_dir, recursive = TRUE); setwd(config$output_dir)
}
cat(red_bold("Working at:", getwd(), "\n"))
void <- file.copy(from = opt$yaml, to = "config.yaml")
register_log <- !interactive() && !grepl("scratch|beegfs", getwd()) && opt$log
log_history <- ""
if(register_log){
  if(opt$verbose) cat(cyan('\n----------------------- Log file\n'))
  log_file <- "_log"; if(opt$verbose) cat("Check log file:", log_file, "\n")
  if(!file.exists(log_file)) file.create(log_file)
  log_file_out <- file(log_file, open = 'wt')
  sink(log_file_out); sink(log_file_out, type = "message")
}
if(opt$verbose) cat('Date and time:\n') ; st.time <- timestamp();

if(opt$verbose){
  cat(cyan('\n----------------------- Parameters\n'))
  str(opt); str(config)
}
if(opt$verbose) cat(cyan('\n----------------------- Getting data\n'))
if(opt$verbose) cat("** Metadata\n")
meta_data = if(!is.null(config$metadata)){
  if(file.exists(config$metadata)) readfile(config$metadata)
}; if(opt$verbose) str(meta_data)

res_file = "scrublet.rds"
if(!file.exists(res_file)){
  # Cell Ranger centric
  tvar <- dir.exists(config$input_matrix[[1]]) && length(config$input_matrix) == 1
  count_file = if(tvar){
    if(grepl("outs|_feature_bc_matrix", config$input_matrix[[1]])){
      config$input_matrix
    }else{
      y <- list.files(path = config$input_matrix, full.names = TRUE)
      paste0(y, "/outs/filtered_feature_bc_matrix")
    }
  }else if(file.exists(config$input_matrix)){
    y <- readfile(config$input_matrix, stringsAsFactors = FALSE)
    z <- apply(X = y, MARGIN = 2, FUN = function(x) sum(file.exists(x)) )
    y[, max(which(z == max(z)))[1]]
  }else{ config$input_matrix }
  tvar <- paste0(config$filters$lib, collapse = "|")
  if(tvar != "") count_file <- count_file[grepl(tvar, count_file)]
  if(!is.null(config$filters$exc))
    count_file <- count_file[!grepl(config$filters$exc, count_file)]
  z <- sapply(config$filters$lib, function(i) which(grepl(i, count_file)) )
  count_file <- count_file[unique(unlist(c(z)))]

  if(opt$verbose) cat("** Matrix data\n")
  # count_data = Seurat::Read10X(count_file) # Using Seurat - might change
  count_data = get_source_data(
    count_file,
    metadata = config$metadata,
    verbose = opt$verbose
  )
  if(!all(rownames(meta_data) %in% colnames(count_file)))
    meta_data <- count_data$mdata
  count_data <- count_data$edata; gc()
  colnames(count_data) = gsub(
    pattern = "([0-9]{1,})_([ACTG]{1,})(.).*",
    replacement = "\\2\\3\\1",
    x = colnames(count_data)
  )
  if(opt$verbose) str(count_data)

  if(length(dim(meta_data)) > 1){
    tvar <- mean(rownames(meta_data) %in% colnames(count_data)) * 100
    cat_col = if(tvar < 1) red_bold else c
    cat(cat_col("Metadata samples found in matrix:", round(tvar, 2), "%\n"))
  }

  library_id = get_library_id(meta_data)
  if(!is.null(library_id)){
    tvar <- table(meta_data[, library_id], useNA = 'always')
    if(requireNamespace("reshape2", quietly = TRUE)) tvar <- reshape2::melt(tvar)
    cat("Libraries\n"); print(tvar)
  }

  if(opt$verbose) cat(red_bold('\n----------------------- Identifying doublets\n'))

  mscrublet = doublets_scrublet(
    edata = count_data,
    mdata = meta_data,
    by = library_id,
    score_threshold = config$filters$score,
    expected_doublet_rate = config$filters$doublet_rate,
    do_parallel = config$do_parallel, #TRUE, #config$do_parallel,
    verbose = opt$verbose
  )
  if(is.null(mscrublet$nCount_RNA)){
    if(opt$verbose) cat("Adding counts column\n")
    csums = if(requireNamespace("Matrix", quietly = TRUE)) Matrix::colSums else base::colSums
    mscrublet$nCount_RNA <- csums(count_data[, rownames(mscrublet)])
  }; tvar <- c("nCount_RNA", "doublet_scores", "predicted_doublets")
  if(opt$verbose) print(summary(mscrublet[, tvar]))
  saveRDS(mscrublet, file = res_file)
}else{
  if(opt$verbose) cat('Existing results\n')
  mscrublet = readRDS(res_file)
}

if(opt$verbose) cat(red_bold('\n----------------------- Plotting results\n'))
if(opt$verbose) str(mscrublet[, !colnames(mscrublet) %in% colnames(meta_data)])
pdf("0_scrublet_threshold.pdf", width = 10)
print(doublets_scrublet_plot(mscrublet, score_threshold = config$filters$score))
graphics.off()

library_id = get_library_id(mscrublet)
if(!is.null(library_id)){
  mscrublet$llibrary = mscrublet[, library_id]
}else{ mscrublet$llibrary = "L1" }
if(length(table(mscrublet$llibrary)) > 1){
  pdf("0_scrublet_threshold_library.pdf", width = 20)
  print(doublets_scrublet_plot(
    mscrublet, by = library_id,
    score_threshold = config$filters$score
  ))
  graphics.off()
}

pass_score <- mscrublet$doublet_scores >= config$filters$score
p0 <- ggplot(
  data = mscrublet,
  mapping = aes_string(x = "doublet_scores", y = "nCount_RNA")
) + geom_point(data = mscrublet[!pass_score, ], size = 0.8, color = "#BEBEBE") +
  geom_point(data = mscrublet[pass_score, ], size = 0.8, color = "black") +
  labs(x = "Doublet score", y = expression(log[10]*"(Count + 1)")) + theme_minimal()
if(requireNamespace("ggpubr", quietly = TRUE))
  p0 <- p0 + ggpubr::stat_cor(method="pearson")
pdf("1_scrublet_vs_counts.pdf")
print(p0)
graphics.off()

if(all(config$vars$class %in% colnames(mscrublet))){
  if(opt$verbose) cat(cyan('\n----------------------- Extra plots\n'))
  source("/home/ciro/scripts/clustering/R/plotting.R")
  source("https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/filters.R")

  tvar = as.character(mscrublet[, config$vars$class[1]])
  tvar[is.na(tvar)] = "Unclassified"
  mscrublet$Method = paste0(
    ifelse(pass_score, "Scrublet-", ""),
    ifelse(grepl(pattern = "Doublet", x = tvar, ignore.case = TRUE), "Hashtag", "")
  )
  mscrublet$Method = gsub("^\\-", "", gsub("\\-$", "", mscrublet$Method))
  mscrublet$Method[mscrublet$Method == ""] = "None"
  table(mscrublet$Method, useNA = "always")
  # cap = max(c(config$filters$score, 0.7))
  # mscrublet$doublet_scores[mscrublet$doublet_scores > cap] <- cap

  p <- simple_violin(mscrublet, config$vars$class[1], "doublet_scores") +
    labs(x = "Hashtag classification", y = "Scrublet score") +
    guides(color = FALSE, fill = FALSE) + theme_minimal() +
    geom_hline(yintercept = config$filters$score, linetype = "dashed", color = "gray50")
  pdf("2_scrublet_hashtag_classification_scores.pdf")
  print(p)
  graphics.off()

  p <- simple_violin(mscrublet[mscrublet$Method != "None", ], "Method", "nCount_RNA") +
    guides(color = FALSE, fill = FALSE) + theme_minimal()
  pdf("2_scrublet_hashtag_classification_counts.pdf")
  print(p)
  graphics.off()

  if(is.null(config$vars$dim))
    config$vars$dim = grep("umap", colnames(mscrublet), value = TRUE, ignore.case = TRUE)
  if(length(config$vars$dim) == 1)
    config$vars$dim = grep(config$vars$dim, colnames(mscrublet), value = TRUE, ignore.case = TRUE)
  if(length(config$vars$dim) > 1){
    if(opt$verbose){
      cat("Dimentionality reduction\n"); str(config$vars$dim)
    }
    p1 <- ggplot(
      data = mscrublet,
      mapping = aes_string(x = config$vars$dim[1], y = config$vars$dim[2])
    ) + geom_point(data = mscrublet[!pass_score, ], size = 0.8, color = "grey") +
      geom_point(data = mscrublet[pass_score, ], size = 0.8, color = "black") +
      labs(x = "Dim 1", y = "Dim 2") + theme_minimal()
    p2 <- ggplot(
      data = mscrublet,
      mapping = aes_string(x = config$vars$dim[1], y = config$vars$dim[2], color = "doublet_scores")
    ) + geom_point(size = 0.8) + theme_minimal() +
      labs(x = "Dim 1", y = "Dim 2", color = "Score")
    if(requireNamespace("viridis", quietly = TRUE))
      p2 <- p2 + viridis::scale_color_viridis(direction=-1)
    pdf("3_scrublet_dim_reduction_scored.pdf", width = 15)
    print(p1 | p2)
    graphics.off()

    p3 <- ggplot(
      data = mscrublet,
      mapping = aes_string(x = config$vars$dim[1], y = config$vars$dim[2], color = "Method")
    ) + geom_point(size = 0.8, color = "#BEBEBE") + theme_minimal() +
      geom_point(data = mscrublet[mscrublet$Method!="None", ], size = 0.8) +
      scale_color_manual(
        values = c(
          None = "#BEBEBE", Hashtag = "#ffa500",
          Scrublet = "#d8db00", "Scrublet-Hashtag" = "#b80704"
        )
      ) + labs(x = "UMAP 1", y = "UMAP 2") +
      guides(colour = guide_legend(override.aes = list(size = 6)))
    pdf("3_scrublet_dim_reduction_methods.pdf", width = 9)
    print(p3)
    graphics.off()
  }
}

if(opt$verbose){
  cat(cyan('\n\n***********************************************************\n'))
  cat('Starting time:\n'); cat(st.time, '\n')
  cat('Finishing time:\n'); timestamp()
  cat('***********************************************************\n')
  cat('SESSION INFO:\n'); print(sessionInfo()); cat("\n")
  cat('Pipeline finished successfully\n')
}

if(register_log){ sink(type = "message"); sink() }
