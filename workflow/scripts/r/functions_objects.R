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
