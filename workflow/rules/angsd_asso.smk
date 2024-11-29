rule create_angsd_asso_ybin_file:
    input:
        samples = config["samples"],
        bams = rules.create_bam_list_allSamples_allSites.output
    output:
        f"{PROGRAM_RESOURCE_DIR}/angsd_habitats.ybin"
    conda: "../envs/picmin.yaml"
    script:
        "../scripts/r/create_angsd_asso_ybin_file.R"

rule angsd_asso_freq:
    """
    Perform association analysis in ANGSD with binary variable (Urban vs. Rural)
    """
    input:
        ybin = rules.create_angsd_asso_ybin_file.output,
        bams = rules.create_bam_list_allSamples_allSites.output,
        ref = rules.copy_ref.output,
        ref_idx = rules.samtools_index_ref.output
    output:
        asso = f'{ANGSD_DIR}/asso/allSamples/{{chrom}}/{{chrom}}_allSamples_asso_freq.lrt0.gz'
    log: f"{LOG_DIR}/angsd_asso/{{chrom}}_angsd_asso_freq.log"
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    params:
        out = f'{ANGSD_DIR}/asso/allSamples/{{chrom}}/{{chrom}}_allSamples_freq'
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 80000,
        runtime = 1440
    shell:
        """
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doMajorMinor 1 \
            -SNP_pval 1e-6 \
            -doMaf 1 \
            -minMaf 0.05 \
            -baq 2 \
            -ref {input.ref} \
            -minQ 20 \
            -minMapQ 30 \
            -doAsso 1 \
            -Pvalue 1 \
            -yBin {input.ybin} \
            -r {wildcards.chrom} \
            -bam {input.bams} 2> {log}
        """

rule angsd_asso_lg:
    """
    Perform association analysis in ANGSD with binary variable (Urban vs. Rural)
    """
    input:
        ybin = rules.create_angsd_asso_ybin_file.output,
        bams = rules.create_bam_list_allSamples_allSites.output,
        ref = rules.copy_ref.output,
        ref_idx = rules.samtools_index_ref.output
    output:
        asso = f'{ANGSD_DIR}/asso/allSamples/{{chrom}}/{{chrom}}_allSamples_asso_lg.lrt0.gz'
    log: f"{LOG_DIR}/angsd_asso/{{chrom}}_angsd_asso_lg.log"
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    params:
        out = f'{ANGSD_DIR}/asso/allSamples/{{chrom}}/{{chrom}}_allSamples_asso_lg'
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 80000,
        runtime = 1440
    shell:
        """
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doMajorMinor 1 \
            -SNP_pval 1e-6 \
            -doMaf 1 \
            -minMaf 0.05 \
            -baq 2 \
            -ref {input.ref} \
            -minQ 20 \
            -minMapQ 30 \
            -doAsso 4 \
            -Pvalue 1 \
            -yBin {input.ybin} \
            -doPost 1 \
            -r {wildcards.chrom} \
            -bam {input.bams} 2> {log}
        """


rule angsd_asso_done:
    input:
        expand(rules.angsd_asso_freq.output, chrom=CHROMOSOMES[0]),
        expand(rules.angsd_asso_lg.output, chrom=CHROMOSOMES[0])
    output:
        f"{ANGSD_DIR}/angsd_asso.done"
    shell:
        """
        touch {output}
        """
