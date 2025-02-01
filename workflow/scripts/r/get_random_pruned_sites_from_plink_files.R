library(tidyverse)

pos <- read_delim(snakemake@input[["pos"]], delim="\t", col_names=FALSE) %>%
    mutate(X1 = str_replace(X1, ":", "_"))
tped <- read_delim(snakemake@input[["tped"]], delim=" ", col_names=FALSE)

set.seed(42)
random_ids <- pos %>%
    sample_n(as.numeric(snakemake@params[["nsites"]])) %>%
    pull(X1)

random_tped <- tped %>%
    filter(X2 %in% random_ids) %>%
    write_delim(., snakemake@output[["tped"]], delim=" ", col_names=FALSE)

read_delim(snakemake@input[["tfam"]], delim=" ", col_names=FALSE) %>%
    write_delim(snakemake@output[["tfam"]], delim=" ", col_names=FALSE)

