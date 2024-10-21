###############
#### SETUP ####
###############

# Load required packages
library(tidyverse)
library(poolr)
library(PicMin)
library(ggnewscale)
source("./scripts/r/functions_objects.R")

#########################
#### PICMIN ANALYSIS ####
#########################

# Analysis here largely follows https://github.com/TBooker/PicMin/blob/main/vignettes/Arabidopsis-vignette.Rmd

all_cities <- read_csv(snakemake@input[["stats"]]) %>%
    filter(Chr == snakemake@wildcards[["chrom"]]) %>%
    dplyr::select(city, winID, mean_c2) %>%
    group_by(city) %>%
    mutate(emp = PicMin:::EmpiricalPs(mean_c2, large_i_small_p = TRUE)) %>%
    dplyr::select(city, winID, emp) 

all_cities_wide <- all_cities %>%
    spread(key = city, value = emp) %>%
    column_to_rownames("winID")

n_lins <- 26  # Total number of cities
n <- 26 

# Run 10,000 replicate simulations of this situation and build the correlation matrix for the order statistics from them
emp_p_null_dat <- t(replicate(40000, PicMin:::GenerateNullData(1.0, n, 0.5, 3, 10000)))

# Calculate the order statistics" p-values for each simulation
emp_p_null_dat_unscaled <- t(apply(emp_p_null_dat, 1, PicMin:::orderStatsPValues))

# Use those p-values to construct the correlation matrix
null_pmax_cor_unscaled <- cor(emp_p_null_dat_unscaled)

# Select the loci that have data for exactly `n_lins` -`n` lineages
lins_p_n <- as.matrix(all_cities_wide[rowSums(is.na(all_cities_wide)) == n_lins - n, ])

# Make some containers for the PicMin results
resulting_p <- rep(-1, nrow(lins_p_n))
resulting_n <- rep(-1, nrow(lins_p_n))

num_reps <- 1000000 # This is an important parameter - the larger the better, but larger values mean longer run times.

# For each of the lines in the dataframe, perform PicMin
for (i in seq_len(nrow(lins_p_n))) {
    test_result <- PicMin:::PicMin(na.omit(lins_p_n[i, ]),
                                   null_pmax_cor_unscaled,
                                   numReps = num_reps)
    # Store the p-value
    resulting_p[i] <- test_result$p
    resulting_n[i] <- test_result$config_est
}

lins_p_n_result <- data.frame(
    numLin = n,
    p = resulting_p,
    q = p.adjust(resulting_p, method = "fdr"),
    n_est = resulting_n,
    locus = row.names(lins_p_n)
)

picmin_results <- cbind(lins_p_n_result,
                        read.csv(
                            text = row.names(lins_p_n),
                            header = FALSE,
                            sep = ":",
                            col.names = c("chrom", "win_center")
                        )) %>%
    mutate(is_outlier = ifelse(-log10(q) >= -log10(snakemake@params[["qval_cut"]]), 1, 0)) %>%
    remap_chr_names()

write_csv(picmin_results, snakemake@output[["picmin"]])

