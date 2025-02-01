library(tidyverse)

ybin <- read_table(snakemake@input[["ybin"]], col_names=FALSE)
fam <- read_delim(snakemake@input[["fam"]], delim=" ", col_names=FALSE) %>%
    mutate(X2 = 1:n()) %>%
    mutate(X6 = ybin$X1)

write_delim(fam, snakemake@output[[1]], delim=" ", col_names = FALSE)
