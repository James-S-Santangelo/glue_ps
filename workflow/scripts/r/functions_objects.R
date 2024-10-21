my_theme <- theme_classic() +
            theme(axis.text = element_text(size = 15),
                  axis.title = element_text(size = 17),
                  legend.title = element_text(size = 15),
                  legend.text = element_text(size = 13)) +
            theme(panel.background = element_rect(fill = "transparent", colour = NA_character_), 
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.background = element_rect(fill = "transparent", colour = NA_character_),
                  legend.background = element_rect(fill = "transparent"),
                  legend.box.background = element_rect(fill = "transparent"),
                  legend.key = element_rect(fill = "transparent"))

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

# Function to load windowed C2 estimates from per-city BayPass runs
load_windowed_c2 <- function(path){

    city <- str_extract(basename(path), "^\\w+(?=_windowed)")
    df <- suppressMessages(read_delim(path, delim="\t")) %>%
        mutate(city = city)
    return(df)
}

remap_chr_names <- function(df){
    df_out <- df %>% 
    mutate(chrom = case_when(chrom == 'Chr01_Occ' ~ 1,
           chrom == 'Chr01_Pall' ~ 2,
           chrom == 'Chr02_Occ' ~ 3,
           chrom == 'Chr02_Pall' ~ 4,
           chrom == 'Chr03_Occ' ~ 5,
           chrom == 'Chr03_Pall' ~ 6,
           chrom == 'Chr04_Occ' ~ 7,
           chrom == 'Chr04_Pall' ~ 8,
           chrom == 'Chr05_Occ' ~ 9,
           chrom == 'Chr05_Pall' ~ 10,
           chrom == 'Chr06_Occ' ~ 11,
           chrom == 'Chr06_Pall' ~ 12,
           chrom == 'Chr07_Occ' ~ 13,
           chrom == 'Chr07_Pall' ~ 14,
           chrom == 'Chr08_Occ' ~ 15,
           chrom == 'Chr08_Pall' ~ 16))
    return(df_out)
}
