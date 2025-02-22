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

rule plink_tped_to_ped:
    input:
        tped = rules.angsd_geno_plink.output.tped,
        tfam = rules.angsd_geno_plink.output.tfam
    output:
        ped = f"{PROGRAM_RESOURCE_DIR}/plink/{{chrom}}_allSamples.ped",
        nosex = f"{PROGRAM_RESOURCE_DIR}/plink/{{chrom}}_allSamples.nosex",
        map = f"{PROGRAM_RESOURCE_DIR}/plink/{{chrom}}_allSamples.map"
    log: f"{LOG_DIR}/plink/tped_to_ped/{{chrom}}.log"
    conda: "../envs/gemma.yaml"
    params:
        out= f"{PROGRAM_RESOURCE_DIR}/plink/{{chrom}}_allSamples",
    shell:
        """
        plink --allow-extra-chr \
            --tfam {input.tfam} \
            --tped {input.tped} \
            --out {params.out} \
            --recode &> {log}
        """

rule plink_generate_bed:
    input:
        ped = rules.plink_tped_to_ped.output.ped,
        map = rules.plink_tped_to_ped.output.map
    output:
        bim = f"{PROGRAM_RESOURCE_DIR}/plink/{{chrom}}_allSamples.bim",
        fam = f"{PROGRAM_RESOURCE_DIR}/plink/{{chrom}}_allSamples_noID_noPheno.fam",
        bed = f"{PROGRAM_RESOURCE_DIR}/plink/{{chrom}}_allSamples.bed"
    log: f"{LOG_DIR}/plink/generate_bed/{{chrom}}.log"
    conda: "../envs/gemma.yaml"
    params:
        out= f"{PROGRAM_RESOURCE_DIR}/plink/{{chrom}}_allSamples",
    shell:
        """
        plink --allow-extra-chr \
            --ped {input.ped} \
            --map {input.map} \
            --out {params.out} \
            --make-bed &> {log} &&

        mv {params.out}.fam {output.fam}
        """

rule add_indID_and_phenotypes:
    input:
        fam = rules.plink_generate_bed.output.fam,
        ybin = rules.create_angsd_asso_ybin_file.output
    output:
        fam = f"{PROGRAM_RESOURCE_DIR}/plink/{{chrom}}_allSamples.fam",
    conda: "../envs/r.yaml"
    script:
        "../scripts/r/add_indID_and_phenotypes.R"

rule get_random_pruned_sites_from_plink_files:
    input:
        pos = rules.prune_graph.output.pos,
        tped = rules.angsd_geno_plink.output.tped,
        tfam = expand(rules.add_indID_and_phenotypes.output.fam, chrom="Chr01_Occ")
    output:
        tped = f"{PROGRAM_RESOURCE_DIR}/plink/random/{{chrom}}_random_pruned.tped",
        tfam = f"{PROGRAM_RESOURCE_DIR}/plink/random/{{chrom}}_random_pruned.tfam"
    conda: "../envs/r.yaml"
    params:
        nsites = 6250
    script:
        "../scripts/r/get_random_pruned_sites_from_plink_files.R"

rule concat_random_tped_files:
    input:
        tped = expand(rules.get_random_pruned_sites_from_plink_files.output.tped, chrom=CHROMOSOMES),
        tfam = expand(rules.get_random_pruned_sites_from_plink_files.output.tfam, chrom="Chr01_Occ")
    output:
        tped = f"{PROGRAM_RESOURCE_DIR}/plink/random/allChroms_random_pruned.tped",
        tfam = f"{PROGRAM_RESOURCE_DIR}/plink/random/allChroms_random_pruned.tfam"
    shell:
        """
        cat {input.tped} > {output.tped}
        mv {input.tfam} {output.tfam}
        """

rule plink_random_tped_to_ped:
    input:
        tped = rules.concat_random_tped_files.output.tped,
        tfam = rules.concat_random_tped_files.output.tfam
    output:
        ped = f"{PROGRAM_RESOURCE_DIR}/plink/random/allChroms_allSamples.ped",
        nosex = f"{PROGRAM_RESOURCE_DIR}/plink/random/allChroms_allSamples.nosex",
        map = f"{PROGRAM_RESOURCE_DIR}/plink/random/allChroms_allSamples.map"
    log: f"{LOG_DIR}/plink/tped_to_ped/allChroms_random.log"
    conda: "../envs/gemma.yaml"
    params:
        out= f"{PROGRAM_RESOURCE_DIR}/plink/random/allChroms_allSamples",
    shell:
        """
        plink --allow-extra-chr \
            --tfam {input.tfam} \
            --tped {input.tped} \
            --out {params.out} \
            --recode &> {log}
        """

rule plink_random_generate_bed:
    input:
        ped = rules.plink_random_tped_to_ped.output.ped,
        map = rules.plink_random_tped_to_ped.output.map
    output:
        bim = f"{PROGRAM_RESOURCE_DIR}/plink/random/allChroms_allSamples.bim",
        fam = f"{PROGRAM_RESOURCE_DIR}/plink/random/allChroms_allSamples.fam",
        bed = f"{PROGRAM_RESOURCE_DIR}/plink/random/allChroms_allSamples.bed"
    log: f"{LOG_DIR}/plink/generate_bed/allChroms_random.log"
    conda: "../envs/gemma.yaml"
    params:
        out= f"{PROGRAM_RESOURCE_DIR}/plink/random/allChroms_allSamples",
    shell:
        """
        plink --allow-extra-chr \
            --ped {input.ped} \
            --map {input.map} \
            --out {params.out} \
            --make-bed &> {log}
        """

rule gemma_estimate_relatedness_matrix:
    input:
        rules.plink_random_generate_bed.output
    output:
        f"{GEMMA_DIR}/relatedness_matrix_random100K.cXX.txt"
    log: f"{LOG_DIR}/gemma/relatedness_matrix.log"
    conda: "../envs/gemma.yaml"
    params:
        outdir = GEMMA_DIR,
        prefix = "relatedness_matrix_random100K",
        bfile = f"{PROGRAM_RESOURCE_DIR}/plink/random/allChroms_allSamples"
    shell:
        """
        gemma -bfile {params.bfile} \
            -outdir {params.outdir} \
            -o {params.prefix} \
            -gk 1 \
            -miss 1.0 \
            -notsnp 2> {log}
        """

rule gemma_ulmm:
    input:
        mat = rules.gemma_estimate_relatedness_matrix.output,
        plink = rules.add_indID_and_phenotypes.output,
    output:
        f"{GEMMA_DIR}/{{chrom}}.assoc.txt"
    log: f"{LOG_DIR}/gemma/{{chrom}}_ulmm.log"
    conda: "../envs/gemma.yaml"
    params:
        outdir = GEMMA_DIR,
        prefix = "{chrom}",
        bfile = f"{PROGRAM_RESOURCE_DIR}/plink/{{chrom}}_allSamples"
    shell:
        """
        gemma -bfile {params.bfile} \
            -outdir {params.outdir} \
            -o {params.prefix} \
            -k {input.mat} \
            -lmm 4 \
            -miss 1.0 \
            -notsnp 2> {log}
        """
        
rule gemma_done:
    input:
        expand(rules.gemma_ulmm.output, chrom=CHROMOSOMES)
    output:
        f"{GEMMA_DIR}/gemma.done"
    shell:
        """
        touch {output}
        """


