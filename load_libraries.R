if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

required_packages = c(
  "Heatplus", "vegan", "viridis", "data.table", "ggplot2", "dplyr", "tidyr",
  "GenomicRanges", "minet", "RColorBrewer", "ggsci", "ComplexHeatmap", "corrplot",
  "tidyverse", "scales", "ggthemes", "Rtsne", "circlize", "GSA", "ggh4x", "ComplexUpset",
  "readxl", "PRROC", "e1071", "randomForest", "stringr", "igraph", "reshape2",
  "WGCNA", "flashClust"
)

missing_packages = required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(missing_packages)) {
  message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
  BiocManager::install(missing_packages, ask = FALSE)
}

invisible(lapply(required_packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))
