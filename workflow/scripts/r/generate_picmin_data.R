###############
#### SETUP ####
###############

# Load required packages
library(tidyverse)
library(PicMin)

# Function to load per-city windowed Fst
load_windowed_fst <- function(path){

    city <- str_extract(basename(path), "^\\w+(?=_Chr)")
    colnames <- c("region", "Chr", "WinCenter", "nSites_fst", "fst")
    df <- suppressMessages(read_delim(path, delim = '\t', skip = 1, col_names = colnames)) %>%
        mutate(fst = ifelse(fst < 0, 0, fst)) %>%
        mutate(city = city)
    return(df)

}

# Function to load per-city windowed thetas
load_windowed_thetas <- function(path){

    city <- str_extract(basename(path), "^\\w+(?=_Chr)")
    habitat <- str_extract(basename(path), "(?<=allSites_)\\w+(?=_win)")
    df <- suppressMessages(read_delim(path, delim = '\t')) %>%
        rename(region = 1,
               "nSites_theta" = "nSites") %>%
        mutate(city = city, habitat = habitat) %>%
        mutate(tp_scaled = tP / nSites_theta) %>%
        mutate(Tajima = ifelse(nSites_theta == 0, NA, Tajima)) %>%
        dplyr::select(region, Chr, WinCenter, tp_scaled, Tajima, nSites_theta, city, habitat)
    return(df)

}

# Load per-city windowed fst
all_windowed_fst <- snakemake@input[["fst"]] %>%
    purrr::map_dfr(., load_windowed_fst, .progress = TRUE) %>%
    mutate(region = str_extract(region, pattern = "(?<=\\()\\d+,\\d+(?=\\)$)")) %>% 
    separate(region, into = c("start", "end"), sep = ",")

# Load per-city windowed thetas
all_windowed_thetas <- snakemake@input[["thetas"]] %>%
    purrr::map_dfr(., load_windowed_thetas, .progress = TRUE) %>%
    mutate(region = str_extract(region, pattern = "(?<=\\()\\d+,\\d+(?=\\)$)")) %>% 
    separate(region, into = c("start", "end"), sep = ",")

all_windowed_thetas_wide <- all_windowed_thetas %>%
    pivot_wider(values_from = c("Tajima", "nSites_theta", "tp_scaled"), names_from = "habitat")

all_windowed_stats <- all_windowed_thetas_wide %>%
    left_join(., all_windowed_fst, by = c("Chr", "start", "end", "WinCenter", "city")) %>%
    mutate(winID = paste0(Chr, ":", WinCenter),
           delta_tp_ur = tp_scaled_urban - tp_scaled_rural,
           delta_td_ur = Tajima_urban - Tajima_rural)

write_csv(all_windowed_stats, snakemake@output[["stats"]])

# Generate subsetted datasets with empirical P-values from Picmin
all_cities_df <- all_windowed_stats %>%
    group_by(city) %>%
    mutate(emp_fst = PicMin:::EmpiricalPs(fst, large_i_small_p = TRUE),
           emp_tp = PicMin:::EmpiricalPs(delta_tp_ur, large_i_small_p = TRUE),
           emp_td = PicMin:::EmpiricalPs(delta_td_ur, large_i_small_p = TRUE)) %>%
    dplyr::select(city, winID, emp_fst, emp_tp, emp_td)

#################################
#### MISSING DATA HISTOGRAMS ####
#################################

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
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))

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
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))

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
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))

ggsave(
    filename = snakemake@output[["hist_td"]],
    plot = missing_data_histogram,
    device = "pdf",
    width = 10,
    height = 10,
    units = "in",
    dpi = 600
)
