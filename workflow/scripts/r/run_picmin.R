###############
#### SETUP ####
###############

# Load required packages
library(tidyverse)
library(poolr)
library(PicMin)
library(ggnewscale)

#########################
#### PICMIN ANALYSIS ####
#########################

# Analysis here largely follows https://github.com/TBooker/PicMin/blob/main/vignettes/Arabidopsis-vignette.Rmd

all_cities_chr1 <- read_csv(snakemake@input[["fst"]]) %>%
    mutate(winID = paste0(Chr, ":", WinCenter)) %>%
    group_by(city) %>%
    mutate(emp = PicMin:::EmpiricalPs(fst, large_i_small_p = TRUE)) %>%
    dplyr::select(city, winID, emp)

all_cities_chr1_wide <- all_cities_chr1 %>%
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

# Select the loci that have data for exactly 7 lineages
lins_p_n <- as.matrix(all_cities_chr1_wide[rowSums(is.na(all_cities_chr1_wide)) == n_lins - n, ])

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
                            col.names = c("scaffold", "win_center")
                        ))

################################
#### PICMIN OUTLIER SUMMARY ####
################################

data_cum <- picmin_results %>%
    group_by(scaffold) %>%
    summarise(max_WinCenter = max(win_center)) %>%
    mutate(WinCenter_add = lag(cumsum(max_WinCenter), default = 0)) %>%
    dplyr::select(scaffold, WinCenter_add)

picmin_results_mod <- picmin_results %>%
    inner_join(data_cum, by = "scaffold") %>%
    mutate(WinCenter_cum = win_center + WinCenter_add)

axis_set <- picmin_results_mod %>%
    group_by(scaffold) %>%
    summarize(center = mean(WinCenter_cum))

outliers <- picmin_results_mod %>%
    filter(-log10(q) >= -log10(0.05))
write_csv(outliers, snakemake@output[[1]])

not_outliers <- picmin_results_mod %>%
    filter(!(-log10(q) >= -log10(0.05)))

manhat_plot <- not_outliers %>%
    mutate(
        chrom_cat = case_when(
            scaffold == "Chr01_Occ" ~ "One",
            scaffold == "Chr01_Pall" ~ "Two",
            scaffold == "Chr02_Occ" ~ "One",
            scaffold == "Chr02_Pall" ~ "Two",
            scaffold == "Chr03_Occ" ~ "One",
            scaffold == "Chr03_Pall" ~ "Two",
            scaffold == "Chr04_Occ" ~ "One",
            scaffold == "Chr04_Pall" ~ "Two",
            scaffold == "Chr05_Occ" ~ "One",
            scaffold == "Chr05_Pall" ~ "Two",
            scaffold == "Chr06_Occ" ~ "One",
            scaffold == "Chr06_Pall" ~ "Two",
            scaffold == "Chr07_Occ" ~ "One",
            scaffold == "Chr07_Pall" ~ "Two",
            scaffold == "Chr08_Occ" ~ "One",
            scaffold == "Chr08_Pall" ~ "Two"
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
        yintercept = -log10(0.05),
        color = "grey40",
        linetype = "dashed"
    ) +
    scale_x_continuous(label = axis_set$scaffold, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 2)) +
    labs(color = "# Cities with evidence of selection",
         fill = "# Cities with evidence of selection") +
    ylab(expression(-log[10] * "(q-value)")) + xlab("") +
    theme_classic() +
    theme(
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(2, "cm"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10, face = "bold"),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(
            size = 14,
            angle = 45,
            hjust = 1,
            vjust = 1
        ),
    ) +
    guides(colour = guide_colourbar(title.position = "top", title.hjust = 0.5))

ggsave(
    filename = snakemake@output[[2]],
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
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))

ggsave(
    filename = snakemake@output[[3]],
    plot = outlier_hist,
    device = "pdf",
    width = 8,
    height = 6,
    units = "in",
    dpi = 600
)
