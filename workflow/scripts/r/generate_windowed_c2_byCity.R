library(tidyverse)

# Load files with chromsome and position order for each split
site_order <- suppressMessages(read_delim(snakemake@input[["site_order"]], 
                                          delim = "\t",
                                          col_names = c("chrom", "pos")))

load_c2_stats <- function(path){
    df <- suppressMessages(read_table(path))
    return(df)
}

c2_df <- snakemake@input[["c2"]] %>% 
    purrr::map_dfr(., load_c2_stats) %>% 
    bind_cols(., site_order) %>% 
    mutate(city = snakemake@wildcards[["city"]]) %>%
    arrange(chrom, pos) %>%
    dplyr::select(city, chrom, pos, C2_std)

calculate_windowed_c2 <- function(df, window_size, step){
    chrom <- df %>% pull(chrom) %>% unique()
    winStarts <- seq(from = 0, to = max(df$pos) + window_size, by = step)
    mat <- matrix(0, nrow = length(winStarts), ncol = 9)
    for(i in 1:length(winStarts)){
        start <- winStarts[i]
        end <- start + step
        df_filt <- df %>% filter(pos >= start & pos < end)
        winID <- i
        winCenter <- start + (step / 2)
        mean <- suppressWarnings(mean(df_filt$C2_std))
        max <- suppressWarnings(max(df_filt$C2_std))
        min <- suppressWarnings(min(df_filt$C2_std))
        n <- nrow(df_filt)
        stats <- c(chrom, winID, start, end, winCenter, mean, max, min, n)
        mat[i, ] <- stats
    }
    stats_df <- as.data.frame(mat)
    names(stats_df) <- c("Chr", "winID", "start", "end", "winCenter", "mean_c2", "max_c2", "min_c2", "n")
    return(stats_df)
}

windowed_c2_df <- c2_df %>%
    group_split(chrom) %>%
    purrr::map_dfr(., calculate_windowed_c2, window_size=10000, step=10000)

write_delim(windowed_c2_df, snakemake@output[["win_c2"]], delim = "\t")
