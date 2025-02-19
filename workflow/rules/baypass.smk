############################################
#### ALLELE COUNTS PER-CITY PER-HABITAT ####
############################################

rule angsd_alleleCounts_byCity_byHabitat:
    input:
        bams = rules.create_bam_list_byHabitat_allSites.output,
        sites = lambda w: [x for x in rules.identify_paralogous_snps.output.sites if w.chrom in x],
        sites_idx = rules.index_filtered_snps.output,
        ref = rules.copy_ref.output,
        ref_idx = rules.samtools_index_ref.output
    output:
        mafs = f'{ANGSD_DIR}/snps/byCity/{{city}}/{{chrom}}/{{city}}_{{habitat}}_{{chrom}}_snps.mafs.gz',
        pos = f'{ANGSD_DIR}/snps/byCity/{{city}}/{{chrom}}/{{city}}_{{habitat}}_{{chrom}}_snps.pos.gz',
        counts = f'{ANGSD_DIR}/snps/byCity/{{city}}/{{chrom}}/{{city}}_{{habitat}}_{{chrom}}_snps.counts.gz',
    log: LOG_DIR + '/angsd_alleleCounts_byCity_byHabitat/{city}/{chrom}_{habitat}_byCity_snps.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    params:
        out = f'{ANGSD_DIR}/snps/byCity/{{city}}/{{chrom}}/{{city}}_{{habitat}}_{{chrom}}_snps'
    threads: 2
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        runtime = lambda wildcards, attempt: attempt * 60
    shell:
        """
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doMajorMinor 4 \
            -baq 2 \
            -ref {input.ref} \
            -doCounts 1 \
            -dumpCounts 3 \
            -doMaf 1 \
            -minQ 20 \
            -minMapQ 30 \
            -sites {input.sites} \
            -anc {input.ref} \
            -r {wildcards.chrom} \
            -bam {input.bams} 2> {log}
        """

rule create_alleleCount_files_byCity:
    input:
        perCity_mafs = expand(rules.angsd_alleleCounts_byCity_byHabitat.output.mafs, city=CITIES, chrom=CHROMOSOMES, habitat=HABITATS),
        perCity_bams = expand(rules.create_bam_list_byHabitat_allSites.output, city=CITIES, habitat=HABITATS),
        global_mafs = expand(rules.angsd_snps_allSamples.output.mafs, chrom=CHROMOSOMES),
        global_sites = expand(rules.identify_paralogous_snps.output.sites, chrom=CHROMOSOMES),
        pos = expand(rules.angsd_alleleCounts_byCity_byHabitat.output.pos, city=CITIES, chrom=CHROMOSOMES, habitat=HABITATS),
        counts = expand(rules.angsd_alleleCounts_byCity_byHabitat.output.counts, city=CITIES, chrom=CHROMOSOMES, habitat=HABITATS),
    output:
        perCity_geno = expand(f"{PROGRAM_RESOURCE_DIR}/baypass/allSites/{{city}}/{{city}}.geno", city=CITIES),
        perCity_cont = expand(f"{PROGRAM_RESOURCE_DIR}/baypass/allSites/{{city}}/{{city}}.cf", city=CITIES),
        perCity_pool = expand(f"{PROGRAM_RESOURCE_DIR}/baypass/allSites/{{city}}/{{city}}.poolsize", city=CITIES),
        site_order = f"{PROGRAM_RESOURCE_DIR}/baypass/allSites/site_order.txt",
        miss = f"{PROGRAM_RESOURCE_DIR}/baypass/allSites/missing.sites"
    conda: "../envs/python.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 40000,
        runtime = lambda wildcards, attempt: attempt * 60 
    params:
        maf_prefix = f"{ANGSD_DIR}/snps",
        out_prefix = f"{PROGRAM_RESOURCE_DIR}/baypass/allSites",
        cities = CITIES,
        habitats = HABITATS,
        chrom = CHROMOSOMES
    notebook:
        "../notebooks/create_alleleCount_files_byCity.py.ipynb"

rule create_baypass_global_input_files:
    input:
        geno = expand(rules.create_alleleCount_files_byCity.output.perCity_geno, city=CITIES),
        cont = expand(rules.create_alleleCount_files_byCity.output.perCity_cont, city=CITIES),
        pool = expand(rules.create_alleleCount_files_byCity.output.perCity_pool, city=CITIES)
    output:
        geno = f"{PROGRAM_RESOURCE_DIR}/baypass/allSites/allSamples.geno",
        cont = f"{PROGRAM_RESOURCE_DIR}/baypass/allSites/allSamples.cont",
        pool = f"{PROGRAM_RESOURCE_DIR}/baypass/allSites/allSamples.pool",
    shell:
        """
        paste {input.geno} | tr '\t' ' ' > {output.geno}
        paste {input.cont} | tr '\t' ' ' > {output.cont}
        paste {input.pool} | tr '\t' ' ' > {output.pool}
        """

rule split_baypass_global_input_files:
    input:
        as_geno = rules.create_baypass_global_input_files.output.geno,
        site_order = rules.create_alleleCount_files_byCity.output.site_order,
        perCity_geno = expand(rules.create_alleleCount_files_byCity.output.perCity_geno, city=CITIES)
    output:
        as_geno = expand(f"{PROGRAM_RESOURCE_DIR}/baypass/split_files/allSamples/splits/allSamples_{{n}}.geno", n=BAYPASS_SPLITS),
        site_order = expand(f"{PROGRAM_RESOURCE_DIR}/baypass/split_files/allSamples/site_order/site_order_{{n}}.txt", n=BAYPASS_SPLITS),
        perCity_geno = expand(f"{PROGRAM_RESOURCE_DIR}/baypass/split_files/byCity/{{city}}/{{city}}_{{n}}.geno", city=CITIES, n=BAYPASS_SPLITS)
    conda: "../envs/python.yaml"
    params:
        out_prefix = f"{PROGRAM_RESOURCE_DIR}/baypass/split_files",
        in_prefix = f"{PROGRAM_RESOURCE_DIR}/baypass/allSites",
        cities = CITIES,
        splits = BAYPASS_SPLITS
    notebook:
        "../notebooks/split_baypass_global_input_files.py.ipynb"

#################
#### BAYPASS ####
#################

rule baypass_coreModel_byCity:
    input:
        unpack(get_baypass_coreModel_byCity_input)
    output:
        log = f"{BAYPASS_DIR}/coreModel_byCity/seed{{k}}/{{city}}/split{{n}}/{{city}}_split{{n}}_baypass.log",
        dic = f"{BAYPASS_DIR}/coreModel_byCity/seed{{k}}/{{city}}/split{{n}}/{{city}}_split{{n}}_DIC.out",
        omega_mat = f"{BAYPASS_DIR}/coreModel_byCity/seed{{k}}/{{city}}/split{{n}}/{{city}}_split{{n}}_mat_omega.out",
        beta_sum = f"{BAYPASS_DIR}/coreModel_byCity/seed{{k}}/{{city}}/split{{n}}/{{city}}_split{{n}}_summary_beta_params.out",
        omega_lda = f"{BAYPASS_DIR}/coreModel_byCity/seed{{k}}/{{city}}/split{{n}}/{{city}}_split{{n}}_summary_lda_omega.out",
        pif_sum = f"{BAYPASS_DIR}/coreModel_byCity/seed{{k}}/{{city}}/split{{n}}/{{city}}_split{{n}}_summary_yij_pij.out",
        pi_xtx = f"{BAYPASS_DIR}/coreModel_byCity/seed{{k}}/{{city}}/split{{n}}/{{city}}_split{{n}}_summary_pi_xtx.out",
        cont_out = f"{BAYPASS_DIR}/coreModel_byCity/seed{{k}}/{{city}}/split{{n}}/{{city}}_split{{n}}_summary_contrast.out",
        bf_out = f"{BAYPASS_DIR}/coreModel_byCity/seed{{k}}/{{city}}/split{{n}}/{{city}}_split{{n}}_summary_betai_reg.out"
    log: f"{LOG_DIR}/baypass/coreModel_byCity/{{city}}/{{city}}_split{{n}}_seed{{k}}.log"
    container: "library://james-s-santangelo/baypass/baypass:2.41"
    threads: 2
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        runtime = lambda wildcards, attempt: attempt * 240
    params:
        out_prefix = f"{BAYPASS_DIR}/coreModel_byCity/seed{{k}}/{{city}}/split{{n}}/{{city}}_split{{n}}"
    shell:
        """
        baypass -gfile {input.geno} \
            -outprefix {params.out_prefix} \
            -seed {wildcards.k} \
            -contrastfile {input.cont} \
            -efile {input.cont} \
            -poolsizefile {input.pool} \
            -nthreads {threads} &> {log}
        """

#############
#### WZA ####
#############

rule create_wza_input_files:
    input:
        urban_mafs = lambda w: expand(rules.angsd_alleleCounts_byCity_byHabitat.output.mafs, chrom=CHROMOSOMES, habitat="urban", city=w.city),
        rural_mafs = lambda w: expand(rules.angsd_alleleCounts_byCity_byHabitat.output.mafs, chrom=CHROMOSOMES, habitat="rural", city=w.city),
        cont_out = lambda w: expand(rules.baypass_coreModel_byCity.output.cont_out, n=BAYPASS_SPLITS, city=w.city, k=[42]),
        site_order = expand(rules.split_baypass_global_input_files.output.site_order, n=BAYPASS_SPLITS)
    output:
        wza_input = f"{PROGRAM_RESOURCE_DIR}/wza/{{city}}_wza_input.txt"
    conda: "../envs/r.yaml"
    params:
        win_size = 20000
    notebook:
        "../notebooks/create_wza_input_files.r.ipynb"

rule run_wza:
    input:
        cor = rules.create_wza_input_files.output.wza_input
    output:
        f"{BAYPASS_DIR}/wza_byCity/{{city}}_wza.txt"
    log: f"{LOG_DIR}/wza/{{city}}_wza.log"
    conda: "../envs/python.yaml"
    shell:
        """
        python3 scripts/python/general_WZA_script.py \
            --correlations {input.cor} \
            --summary_stat C2_pval \
            --window win \
            --MAF mean_maf \
            --output {output} \
            --sep "\t" \
            --verbose 2> {log}
        """

##################
#### ANALYSES ####
##################

# rule fmd_and_omega_mat_pca:
#     input:
#         omega_mat = rules.baypass_coreModel_random100K.output.omega_mat
#     output:
#         pca = f"{ANALYSIS_DIR}/baypass/figures/pca_by_continent_and_habitat.pdf",
#     conda: "../envs/r.yaml"
#     params:
#         cities = CITIES,
#         habitats = HABITATS,
#     notebook:
#         "../notebooks/fmd_and_omega_mat_pca.r.ipynb"

rule baypass_outlier_test:
    input:
        cont_out = expand(rules.baypass_coreModel_byCity.output.cont_out, n=BAYPASS_SPLITS, city=CITIES, k=[42]),
        pi_xtx = expand(rules.baypass_coreModel_byCity.output.pi_xtx, n = BAYPASS_SPLITS, city=CITIES, k=[42]),
        site_order = expand(rules.split_baypass_global_input_files.output.site_order, n=BAYPASS_SPLITS),
    output:
        "test.txt",
        c2_outliers = f"{ANALYSIS_DIR}/baypass/baypass_c2_outliers.txt",
        c2_pval_hist = f"{ANALYSIS_DIR}/baypass/figures/baypass_c2_pval_hist.pdf",
        c2_manhat_pdf = f"{ANALYSIS_DIR}/baypass/figures/baypass_c2_manhattan.pdf",
        c2_manhat_png = f"{ANALYSIS_DIR}/baypass/figures/baypass_c2_manhattan.png"
    params:
        qval_cut = 0.05
    conda: "../envs/r.yaml"
    notebook:
        "../notebooks/baypass_outlier_test.r.ipynb"

##############
#### POST ####
##############

rule baypass_done:
    input:
        # expand(rules.generate_windowed_c2_byCity.output, city=CITIES),
        # rules.fmd_and_omega_mat_pca.output,
        # rules.baypass_outlier_test.output,
        expand(rules.run_wza.output, city=CITIES),
        expand(rules.baypass_coreModel_byCity.output, city=CITIES, n=BAYPASS_SPLITS, k=[42]),
    output:
        f"{BAYPASS_DIR}/baypass.done"
    shell:
        """
        touch {output}
        """

