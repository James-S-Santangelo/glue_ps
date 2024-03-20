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

picmin_results <- read_csv(snakemake@input[["results"]])
picmin_results$pooled_q <- p.adjust(picmin_results$p, method = "fdr")
write_csv(picmin_results, snakemake@output[["picmin"]])

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
    filter(-log10(pooled_q) >= -log10(0.05))
write_csv(outliers, snakemake@output[["outlier"]])

not_outliers <- picmin_results_mod %>%
    filter(!(-log10(pooled_q) >= -log10(0.05)))

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
    ggplot(aes(x = WinCenter_cum, y = -log10(pooled_q))) +
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
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))

ggsave(
    filename = snakemake@output[["hist"]],
    plot = outlier_hist,
    device = "pdf",
    width = 8,
    height = 6,
    units = "in",
    dpi = 600
)
