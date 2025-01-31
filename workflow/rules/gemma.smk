rule angsd_geno_plink:
    """
    Call genotypes and export to plink format
    """
    input:
        bams = rules.create_bam_list_allSamples_allSites.output,
        ref = rules.copy_ref.output,
        ref_idx = rules.samtools_index_ref.output,
        sites = lambda w: [x for x in rules.identify_paralogous_snps.output.sites if w.chrom in x],
        sites_idx = rules.index_filtered_snps.output
    output:
        geno = f'{ANGSD_DIR}/geno/allSamples/{{chrom}}_allSamples.geno.gz',
        tped = f'{ANGSD_DIR}/geno/allSamples/{{chrom}}_allSamples.tped',
        tfam = f'{ANGSD_DIR}/geno/allSamples/{{chrom}}_allSamples.tfam'
    log: f"{LOG_DIR}/angsd_geno_plink/{{chrom}}_angsd_geno_plink.log"
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    params:
        out = f'{ANGSD_DIR}/geno/allSamples/{{chrom}}_allSamples'
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        runtime = 1440
    shell:
        """
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doMajorMinor 4 \
            -doMaf 1 \
            -baq 2 \
            -ref {input.ref} \
            -minQ 20 \
            -minMapQ 30 \
            -remove_bads 1 \
            -skipTriallelic 1 \
            -uniqueOnly 1 \
            -only_proper_pairs 1 \
            -doPost 1 \
            -postCutoff 0.95 \
            -doGeno 5 \
            -doPlink 2 \
            -r {wildcards.chrom} \
            -sites {input.sites} \
            -bam {input.bams} 2> {log}
        """

rule gemma_done:
    input:
        expand(rules.angsd_geno_plink.output, chrom=CHROMOSOMES)
    output:
        f"{GEMMA_DIR}/gemma.done"
    shell:
        """
        touch {output}
        """


