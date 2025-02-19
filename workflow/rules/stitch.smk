################################
#### PARAMETER OPTIMIZATION ####
################################

# Run STITCH on single chromosome with different parameter combinations
# Evaluate accuracy and performance and select optimat parameters for genome-wide run

rule get_highCov_reference_sites:
    """
    Generate text file with list of sites in high coverage reference VCF
    """
    input:
        HIGH_COV_VCF
    output:
        f"{PROGRAM_RESOURCE_DIR}/stitch/optimization/Chr01_Occ_highCov_samples.sites"
    shell:
        """
        zgrep -v '^#' {input} | cut -f1,2 > {output}
        """

rule get_overlapping_samples:
    """
    Create text file with GLUE samples from Toronto (which have high cov data as well)
    """
    input:
        config["samples"]
    output:
        f"{PROGRAM_RESOURCE_DIR}/stitch/optimization/toronto_samples.txt"
    shell:
        """
        grep 'Toronto' {input} | cut -f6 > {output}
        """


rule convert_angsd_mafs_to_sites:
    """
    Extract chromosome, position, ref allele, alt allele from ANGSD output for chromosome 1
    """
    input:
        expand(rules.angsd_snps_allSamples.output.mafs, chrom="Chr01_Occ")
    output:
        f"{PROGRAM_RESOURCE_DIR}/stitch/optimization/Chr01_Occ_angsd_snps_allSamples.sites"
    shell:
        """
        zcat {input} | cut -f1,2,3,4 > {output}
        """

rule filter_overlapping_sites:
    """
    Fileter ANGSD SNPs to only those with genotype calls in the high coverence reference set
    """
    input:
        angsd_sites = rules.convert_angsd_mafs_to_sites.output,
        high_cov_sites = rules.get_highCov_reference_sites.output
    output:
        f"{PROGRAM_RESOURCE_DIR}/stitch/optimization/Chr01_Occ_stitch.sites"
    conda: "../envs/python.yaml"
    script:
        "../scripts/python/filter_overlapping_sites.py"

rule filter_subset_vcf:
    """
    Removes samples from high coverage reference VCF that were not sequenced as part of GLUE.
    Filters to only include sites overlapping ANGSD SNPs
    """
    input:
        vcf = HIGH_COV_VCF,
        samples = rules.get_overlapping_samples.output,
        sites = rules.filter_overlapping_sites.output
    output:
        f"{PROGRAM_RESOURCE_DIR}/stitch/optimization/Chr01_Occ_subset_filtered.vcf.gz"
    log: f"{LOG_DIR}/keep_overlapping_samples/keep_overlapping_samples.log"
    conda: "../envs/stitch.yaml"
    shell:
        """
        bcftools view -O z -o {output} \
            --samples-file {input.samples} \
            --force-samples \
            --targets-file <(cut -f1,2 {input.sites}) \
            {input.vcf} 2> {log}
        """

rule generate_dosage_matrix:
    """
    Convert high coverage reference VCF to dosage matrix for STITCH
    """
    input:
        vcf = rules.filter_subset_vcf.output
    output:
        f"{PROGRAM_RESOURCE_DIR}/stitch/optimization/Chr01_Occ_dosages.mat"
    container: "docker://ghcr.io/thewanglab/algatr"
    script:
        "../scripts/r/generate_dosage_matrix.R"

rule stitch_done:
    input:
        rules.filter_overlapping_sites.output,
        rules.generate_dosage_matrix.output
    output:
        f"{STITCH_DIR}/stitch.done"
    shell:
        """
        touch {output}
        """
