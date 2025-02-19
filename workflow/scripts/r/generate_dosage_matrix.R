library(algatr)
library(readr)

vcf <- vcfR::read.vcfR(snakemake@input[["vcf"]])
dos_mat <- vcf_to_dosage(vcf)
print(snakemake@output[[1]])
write_delim(as.data.frame(t(dos_mat)), snakemake@output[[1]], delim="\t")
