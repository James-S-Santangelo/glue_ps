###############
#### SETUP ####
###############

# Load required packages
library(tidyverse)
library(PicMin)
source("./scripts/r/functions_objects.R")

# Load per-city windowed fst
all_windowed_fst <- snakemake@input[["fst"]] %>%
    purrr::map_dfr(., load_windowed_fst, .progress = TRUE) %>%
    mutate(region = str_extract(region, pattern = "(?<=\\()\\d+,\\d+(?=\\)$)")) %>% 
    separate(region, into = c("start", "end"), sep = ",") %>%
    mutate(start = as.character(start), end = as.character(end))

# Load all windowed C2 estimates
all_windowed_c2 <- snakemake@input[["c2"]] %>%
    purrr::map_dfr(., load_windowed_c2, .progress = TRUE) %>%
    rename("WinCenter" = "winCenter") %>%
    mutate(start = as.character(start), end = as.character(end))

# Load per-city windowed thetas
all_windowed_thetas <- snakemake@input[["thetas"]] %>%
    purrr::map_dfr(., load_windowed_thetas, .progress = TRUE) %>%
    mutate(region = str_extract(region, pattern = "(?<=\\()\\d+,\\d+(?=\\)$)")) %>% 
    separate(region, into = c("start", "end"), sep = ",") %>%
    mutate(start = as.character(start), end = as.character(end))

all_windowed_thetas_wide <- all_windowed_thetas %>%
    pivot_wider(values_from = c("Tajima", "nSites_theta", "tp_scaled"), names_from = "habitat")

all_windowed_stats <- all_windowed_thetas_wide %>%
    left_join(., all_windowed_fst, by = c("Chr", "start", "end", "WinCenter", "city")) %>%
    left_join(., all_windowed_c2, by = c("Chr", "start", "end", "WinCenter", "city")) %>%
    mutate(winID = paste0(Chr, ":", WinCenter),
           delta_tp_ur = tp_scaled_urban - tp_scaled_rural,
           delta_td_ur = Tajima_urban - Tajima_rural,
           abs_delta_tp_ur = abs(delta_tp_ur),
           abs_delta_td_ur = abs(delta_td_ur))

write_csv(all_windowed_stats, snakemake@output[["stats"]])

#################################
#### MISSING DATA HISTOGRAMS ####
#################################

# Generate subsetted datasets with empirical P-values from Picmin
all_cities_df <- all_windowed_stats %>%
    group_by(city) %>%
    mutate(emp_fst = PicMin:::EmpiricalPs(fst, large_i_small_p = TRUE),
           emp_c2 = PicMin:::EmpiricalPs(mean_c2, large_i_small_p = TRUE),
           emp_tp = PicMin:::EmpiricalPs(abs_delta_tp_ur, large_i_small_p = TRUE),
           emp_td = PicMin:::EmpiricalPs(abs_delta_td_ur, large_i_small_p = TRUE)) %>%
    dplyr::select(city, winID, emp_fst, emp_c2, emp_tp, emp_td)

all_cities_wide <- all_cities_df %>%
    dplyr::select(city, winID, emp_fst) %>%
    spread(key = city, value = emp_fst) %>%
    column_to_rownames("winID")

column_names <- colnames(all_cities_wide)
missing_data_histogram <- all_cities_wide %>%
    rowwise() %>%
    mutate(na_count = sum(is.na(c_across(all_of(column_names))))) %>%
    ggplot(aes(x = na_count)) +
    geom_histogram(bins = 26, color = "black", fill = "white") +
    ylab("Number of genomic windows") + xlab("Number of cities with missing data") +
    theme_classic() +
    my_theme

ggsave(
    filename = snakemake@output[["hist_fst"]],
    plot = missing_data_histogram,
    device = "pdf",
    width = 10,
    height = 10,
    units = "in",
    dpi = 600
)

all_cities_wide <- all_cities_df %>%
    dplyr::select(city, winID, emp_c2) %>%
    spread(key = city, value = emp_c2) %>%
    column_to_rownames("winID")

column_names <- colnames(all_cities_wide)
missing_data_histogram <- all_cities_wide %>%
    rowwise() %>%
    mutate(na_count = sum(is.na(c_across(all_of(column_names))))) %>%
    ggplot(aes(x = na_count)) +
    geom_histogram(bins = 26, color = "black", fill = "white") +
    ylab("Number of genomic windows") + xlab("Number of cities with missing data") +
    theme_classic() +
    my_theme

ggsave(
    filename = snakemake@output[["hist_c2"]],
    plot = missing_data_histogram,
    device = "pdf",
    width = 10,
    height = 10,
    units = "in",
    dpi = 600
)

all_cities_wide <- all_cities_df %>%
    dplyr::select(city, winID, emp_tp) %>%
    spread(key = city, value = emp_tp) %>%
    column_to_rownames("winID")

column_names <- colnames(all_cities_wide)
missing_data_histogram <- all_cities_wide %>%
    rowwise() %>%
    mutate(na_count = sum(is.na(c_across(all_of(column_names))))) %>%
    ggplot(aes(x = na_count)) +
    geom_histogram(bins = 26, color = "black", fill = "white") +
    ylab("Number of genomic windows") + xlab("Number of cities with missing data") +
    theme_classic() +
    my_theme

ggsave(
    filename = snakemake@output[["hist_tp"]],
    plot = missing_data_histogram,
    device = "pdf",
    width = 10,
    height = 10,
    units = "in",
    dpi = 600
)

all_cities_wide <- all_cities_df %>%
    dplyr::select(city, winID, emp_td) %>%
    spread(key = city, value = emp_td) %>%
    column_to_rownames("winID")

column_names <- colnames(all_cities_wide)
missing_data_histogram <- all_cities_wide %>%
    rowwise() %>%
    mutate(na_count = sum(is.na(c_across(all_of(column_names))))) %>%
    ggplot(aes(x = na_count)) +
    geom_histogram(bins = 26, color = "black", fill = "white") +
    ylab("Number of genomic windows") + xlab("Number of cities with missing data") +
    theme_classic() +
    my_theme

ggsave(
    filename = snakemake@output[["hist_td"]],
    plot = missing_data_histogram,
    device = "pdf",
    width = 10,
    height = 10,
    units = "in",
    dpi = 600
)
