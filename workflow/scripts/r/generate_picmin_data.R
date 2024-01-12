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

# Load per-city windowed fst, concatenate, and write
windowed_fst_list <- snakemake@input[["fst"]]
all_windowed_fst <- windowed_fst_list %>%
    purrr::map_dfr(., load_windowed_fst) %>%
    mutate(region = str_extract(region, pattern = "(?<=\\()\\d+,\\d+(?=\\)$)")) %>% 
    separate(region, into = c("start", "end"), sep = ",")

write_csv(all_windowed_fst, snakemake@output[[1]])

# Generate subsetted datasets with empirical P-values from Picmin
all_cities_chr1 <- all_windowed_fst %>%
    mutate(winID = paste0(Chr, ":", WinCenter)) %>%
    group_by(city) %>%
    mutate(emp = PicMin:::EmpiricalPs(fst, large_i_small_p = TRUE)) %>%
    dplyr::select(city, winID, emp)

################################
#### MISSING DATA HISTOGRAM ####
################################

all_cities_chr1_wide <- all_cities_chr1 %>%
    spread(key = city, value = emp) %>%
    column_to_rownames("winID")

column_names <- colnames(all_cities_chr1_wide)
missing_data_histogram <- all_cities_chr1_wide %>%
    rowwise() %>%
    mutate(na_count = sum(is.na(c_across(all_of(column_names))))) %>%
    ggplot(aes(x = na_count)) +
    geom_histogram(bins = 26, color = "black", fill = "white") +
    ylab("Number of genomic windows") + xlab("Number of cities with missing data") +
    theme_classic() +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))

ggsave(
    filename = snakemake@output[[2]],
    plot = missing_data_histogram,
    device = "pdf",
    width = 10,
    height = 10,
    units = "in",
    dpi = 600
)
