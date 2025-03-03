library(STITCH)

K <- snakemake@wildcards[["K"]]
nGen <- snakemake@wildcards[["nGen"]]

outname <- paste0("Chr01_Occ_K", K, "_nGen", nGen, ".vcf.gz") 

log <- file(snakemake@log[[1]])
sink(log, append=TRUE)
sink(log, append=TRUE, type="message")

STITCH(
    chr = "Chr01_Occ",
    bamlist = snakemake@input[["bams"]],
    posfile = snakemake@input[["posfile"]],
    genfile = snakemake@input[["dosage"]],
    outputdir = snakemake@params[["outdir"]],
    output_filename = outname,
    nCores = snakemake@threads,
    K = as.integer(K),
    nGen = as.integer(nGen),
    regionStart = 10000000,
    regionEnd = 20000000,
    buffer = 100000,
    method = "diploid",
    switchModelIteration = 59,
    output_format = "bgvcf",
    niterations = 60,
    shuffleHaplotypeIterations = c(5, 10, 15, 20, 25, 30, 35, 40),
    expRate = 5.0,
    plot_shuffle_haplotype_attempts = TRUE,
    shuffle_bin_radius = 100 
)
