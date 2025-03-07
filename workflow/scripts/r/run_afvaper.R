library(afvaper)
library(tidyverse)

# Load frequency matrix
af_mat <- read_delim(snakemake@input[["afs"]], delim = "\t") %>%
    as.data.frame()

# Build vector list
vector_list <- list(c("Albuquerque_rural","Albuquerque_urban"),
                    c("Antwerp_rural","Antwerp_urban"),
                    c("Armidale_rural","Armidale_urban"),
                    c("Athens_rural","Athens_urban"),
                    c("Bogota_rural","Bogota_urban"),
                    c("Buenos_Aires_rural","Buenos_Aires_urban"),
                    c("Calgary_rural","Calgary_urban"),
                    c("Canberra_rural","Canberra_urban"),
                    c("Cape_Town_rural","Cape_Town_urban"),
                    c("Christchurch_rural","Christchurch_urban"),
                    c("Kunming_rural","Kunming_urban"),
                    c("Landshut_rural","Landshut_urban"),
                    c("Linkoping_rural","Linkoping_urban"),
                    c("Loja_rural","Loja_urban"),
                    c("Medellin_rural","Medellin_urban"),
                    c("Memphis_rural","Memphis_urban"),
                    c("Munich_rural","Munich_urban"),
                    c("Palmerston_North_rural","Palmerston_North_urban"),
                    c("Punta_Arenas_rural","Punta_Arenas_urban"),
                    c("Quito_rural","Quito_urban"),
                    c("Sapporo_rural","Sapporo_urban"),
                    c("Tehran_rural","Tehran_urban"),
                    c("Thessaloniki_rural","Thessaloniki_urban"),
                    c("Toronto_rural","Toronto_urban"),
                    c("Vancouver_rural","Vancouver_urban"),
                    c("Warsaw_rural","Warsaw_urban"))

# Name vectors
names(vector_list) <- c("Albuquerque","Antwerp","Armidale","Athens","Bogota",
                        "Buenos_Aires","Calgary","Canberra","Cape_Town","Christchurch",
                        "Kunming","Landshut","Linkoping","Loja","Medellin",
                        "Memphis","Munich","Palmerston_North","Punta_Arenas","Quito",
                        "Sapporo","Tehran","Thessaloniki","Toronto","Vancouver",
                        "Warsaw")

# Determine number of permutations per chromosome
fai <- read_delim(snakemake@input[["fai"]], delim="\t", col_names = FALSE) %>%
    dplyr::select(1,2)
names(fai) <- c("chr", "size")
fai <- fai %>%
   dplyr::filter(!(chr %in% c("Mitochondria", "Plastid")))

# How many permutations do we want in total?
total_perms <- 1000000

# Fetch proportional size of all chromosomes
chr_props <- fai$size/sum(fai$size)
chr_perms <- data.frame(chr=fai$chr,
                        perms=round(chr_props * total_perms))

# Set our window size
window_snps <- 200

# Calculate Allele Frequency Change Vector Matrices
AF_input <- calc_AF_vectors(vcf = af_mat,
                            window_size = window_snps,
                            vectors = vector_list,
                            n_cores = 1,
                            data_type = "freq")

# Get number of permutations for current chromosome
actual_perms <- chr_perms %>%
    filter(chr == snakemake@wildcards[["chrom"]]) %>%
    pull(perms)

# Calculate Allele Frequency Change Vector Matrices
null_input <- calc_AF_vectors(vcf = af_mat,
                              window_size = window_snps,
                              vectors = vector_list,
                              n_cores = 1,
                              null_perms = actual_perms,
                              data_type = "freq")

saveRDS(list(AF_input, null_input), snakemake@output[[1]])
