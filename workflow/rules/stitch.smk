################################
#### PARAMETER OPTIMIZATION ####
################################

# Run STITCH on single chromosome with different parameter combinations
# Evaluate accuracy and performance and select optimal parameters for genome-wide run

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

rule filter_subset_vcf:
    """
    Removes samples from high coverage reference VCF that were not sequenced as part of GLUE.
    Filters to only include sites overlapping ANGSD SNPs
    """
    input:
        vcf = HIGH_COV_VCF,
        samples = rules.get_overlapping_samples.output
    output:
        vcf = f"{PROGRAM_RESOURCE_DIR}/stitch/optimization/Chr01_Occ_subset_filtered.vcf.gz",
        tbi = f"{PROGRAM_RESOURCE_DIR}/stitch/optimization/Chr01_Occ_subset_filtered.vcf.gz.tbi"
    log: f"{LOG_DIR}/keep_overlapping_samples/keep_overlapping_samples.log"
    conda: "../envs/stitch.yaml"
    shell:
        """
        ( bcftools view \
            --samples-file {input.samples} \
            --force-samples  {input.vcf} |
        bcftools filter -i 'F_MISSING <= 0.7' |
            bcftools filter -i 'QUAL >= 30' |
            bcftools filter -i 'AB >= 0.25 & AB <= 0.75 | AB <= 0.01' |
            bcftools filter -i 'SAF > 0 & SAR > 0' |
            bcftools filter -i 'MQM >=30 & MQMR >= 30' |
            bcftools filter -i '((PAIRED > 0.05) & (PAIREDR > 0.05) & (PAIREDR / PAIRED < 1.75 ) & (PAIREDR / PAIRED > 0.25)) | ((PAIRED < 0.05) & (PAIREDR < 0.05))' |
            bcftools filter -i '((AF > 0) & (AF < 1))' -O z -o {output.vcf} &&
        sleep 5
        tabix {output.vcf} ) 2> {log}
        """

rule get_highCov_reference_sites:
    """
    Generate text file with list of sites in high coverage reference VCF
    """
    input:
        rules.filter_subset_vcf.output.vcf
    output:
        f"{PROGRAM_RESOURCE_DIR}/stitch/optimization/Chr01_Occ_highCov_samples.sites"
    shell:
        """
        zgrep -v '^#' {input} | cut -f1,2 > {output}
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

rule get_overlapping_sites_fromVCF:
    """
    Filter high coverage VCF to only those sites that overlap with ANGSD SNPs from GLUE
    """
    input:
        vcf = rules.filter_subset_vcf.output.vcf,
        sites = rules.filter_overlapping_sites.output
    output:
        f"{PROGRAM_RESOURCE_DIR}/stitch/optimization/Chr01_Occ_subset_filtered_overlappingSites.vcf.gz"
    log: f"{LOG_DIR}/get_overlapping_sites_fromVCF/get_overlapping_sites_fromVCF.log"
    conda: "../envs/stitch.yaml"
    shell:
        """
        vcfuniq {input.vcf} | bcftools view -O z -o {output} \
            --targets-file <(cut -f1,2 {input.sites}) 2> {log}
        """

rule generate_dosage_matrix:
    """
    Convert high coverage reference VCF to dosage matrix for STITCH
    """
    input:
        vcf = rules.get_overlapping_sites_fromVCF.output
    output:
        f"{PROGRAM_RESOURCE_DIR}/stitch/optimization/Chr01_Occ_dosages.mat"
    container: "docker://ghcr.io/thewanglab/algatr"
    script:
        "../scripts/r/generate_dosage_matrix.R"

rule run_stitch:
    """
    Run STITCH to generate imputed genotypes
    """
    input:
        bams = rules.create_bam_list_allSamples_allSites.output,
        posfile = rules.filter_overlapping_sites.output,
        dosage = rules.generate_dosage_matrix.output
    output:
        f"{STITCH_DIR}/optimization/K{{K}}_nGen{{nGen}}/Chr01_Occ_K{{K}}_nGen{{nGen}}.vcf.gz"
    log: f"{LOG_DIR}/stitch/optimization/K{{K}}_{{nGen}}.log"
    conda: "../envs/stitch.yaml"
    threads: 6
    params:
        outdir = f"{STITCH_DIR}/optimization/K{{K}}_nGen{{nGen}}"
    script:
        "../scripts/r/run_stitch.R"

rule stitch_done:
    input:
        expand(rules.run_stitch.output, K=[10, 50, 100], nGen=[1000, 10000, 100000])
    output:
        f"{STITCH_DIR}/stitch.done"
    shell:
        """
        touch {output}
        """
