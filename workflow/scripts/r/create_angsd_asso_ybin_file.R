library(tidyverse)

samples <- read_delim(snakemake@input[["samples"]], delim="\t") %>%
    mutate(ybin = ifelse(site == "urban", 1, 0)) %>%
    dplyr::select(sample, ybin)

bams <- read_table(snakemake@input[["bams"]], col_names = c("bam")) %>%
    separate(bam, sep = "/", into = c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10")) %>%
    mutate(sample = str_extract(V10, ".*(?=_merged)")) %>%
    dplyr::select(sample)

ybin_file <- bams %>%
    left_join(samples, by = "sample") %>%
    dplyr::select(ybin)

write_delim(ybin_file, snakemake@output[[1]], col_names = FALSE)
