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

if (snakemake@wildcards[["stat"]] == "fst") {
    all_cities <- read_csv(snakemake@input[["stats"]]) %>%
        dplyr::select(city, winID, fst) %>%
        group_by(city) %>%
        mutate(emp = PicMin:::EmpiricalPs(fst, large_i_small_p = TRUE)) %>%
        dplyr::select(city, winID, emp)
} else if (snakemake@wildcards[["stat"]] == "tp") {
   all_cities <- read_csv(snakemake@input[["stats"]]) %>%
       dplyr::select(city, winID, abs_delta_tp_ur) %>%
       group_by(city) %>%
       mutate(emp = PicMin:::EmpiricalPs(abs_delta_tp_ur, large_i_small_p = TRUE)) %>%
       dplyr::select(city, winID, emp)
} else {
    all_cities <- read_csv(snakemake@input[["stats"]]) %>%
       dplyr::select(city, winID, abs_delta_td_ur) %>%
        group_by(city) %>%
        mutate(emp = PicMin:::EmpiricalPs(abs_delta_td_ur, large_i_small_p = TRUE)) %>%
        dplyr::select(city, winID, emp)
}


all_cities_wide <- all_cities %>%
    spread(key = city, value = emp) %>%
    column_to_rownames("winID")

n_lins <- 26  # Total number of cities
n <- as.numeric(snakemake@wildcards[["n"]])  # Number of cities with data (i.e., not missing)

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

################################
#### PICMIN OUTLIER SUMMARY ####
################################

data_cum <- picmin_results %>%
    group_by(chrom) %>%
    summarise(max_WinCenter = max(win_center)) %>%
    mutate(WinCenter_add = lag(cumsum(max_WinCenter), default = 0)) %>%
    dplyr::select(chrom, WinCenter_add)

picmin_results_mod <- picmin_results %>%
    inner_join(data_cum, by = "chrom") %>%
    mutate(WinCenter_cum = win_center + WinCenter_add)

axis_set <- picmin_results_mod %>%
    group_by(chrom) %>%
    summarize(center = mean(WinCenter_cum))

outliers <- picmin_results_mod %>%
    filter(is_outlier == 1)

not_outliers <- picmin_results_mod %>%
    filter(is_outlier == 0)

manhat_plot <- not_outliers %>%
    mutate(
        chrom_cat = case_when(
            chrom == 1 ~ "One",
            chrom == 2 ~ "Two",
            chrom == 3 ~ "One",
            chrom == 4 ~ "Two",
            chrom == 5 ~ "One",
            chrom == 6 ~ "Two",
            chrom == 7 ~ "One",
            chrom == 8 ~ "Two",
            chrom == 9 ~ "One",
            chrom == 10 ~ "Two",
            chrom == 11 ~ "One",
            chrom == 12 ~ "Two",
            chrom == 13 ~ "One",
            chrom == 14 ~ "Two",
            chrom == 15 ~ "One",
            chrom == 16 ~ "Two"
        )
    ) %>%
    ggplot(aes(x = WinCenter_cum, y = -log10(q))) +
    geom_point(
        shape = 21,
        size = 3,
        aes(fill = chrom_cat, color = chrom_cat),
        show.legend = FALSE
    ) +
    scale_fill_manual(values = c("black", "grey40")) +
    scale_color_manual(values = c("black", "grey40")) +
    new_scale_color() +
    new_scale_fill() +
    geom_point(
        data = outliers,
        shape = 21,
        size = 3,
        aes(color = n_est, fill = n_est)
    ) +
    scale_colour_viridis_c(option = "plasma", direction = -1) +
    scale_fill_viridis_c(option = "plasma", direction = -1) +
    geom_hline(
        yintercept = -log10(snakemake@params[["qval_cut"]]),
        color = "grey40",
        linetype = "dashed"
    ) +
    scale_x_continuous(label = axis_set$chrom, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 2)) +
    labs(color = "# Cities with evidence of selection",
         fill = "# Cities with evidence of selection") +
    ylab(expression(-log[10] * "(q-value)")) + xlab("") +
    theme_classic() +
    my_theme +
    theme(
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(2, "cm"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold")
    ) +
    guides(colour = guide_colourbar(title.position = "top", title.hjust = 0.5))

ggsave(
    filename = snakemake@output[["manhat"]],
    plot = manhat_plot,
    device = "pdf",
    width = 20,
    height = 7,
    units = "in",
    dpi = 600
)

outlier_hist <- ggplot(outliers, aes(x = n_est)) +
    geom_histogram(bins = 21,
                   color = "black",
                   fill = "white") +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 10, 2)) +
    xlab("# of cities with evidence of selection") +
    ylab("# of genomic windows") +
    theme_classic()  +
    my_theme

ggsave(
    filename = snakemake@output[["hist"]],
    plot = outlier_hist,
    device = "pdf",
    width = 8,
    height = 6,
    units = "in",
    dpi = 600
)
