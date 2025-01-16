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
    conda: "../envs/baypass.yaml"
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

#################
#### BAYPASS ####
#################

##################
#### ANALYSES ####
##################

# rule fmd_and_omega_mat_pca:
#     input:
#         omega_mat = rules.baypass_coreModel_random100K.output.omega_mat
#     output:
#         pca = f"{ANALYSIS_DIR}/baypass/figures/pca_by_continent_and_habitat.pdf",
#     conda: "../envs/baypass.yaml"
#     params:
#         cities = CITIES,
#         habitats = HABITATS,
#     notebook:
#         "../notebooks/fmd_and_omega_mat_pca.r.ipynb"

# rule baypass_outlier_test:
#     input:
#         cont_out = expand(rules.baypass_coreModel_allSamples.output.cont_out, n=BAYPASS_SPLITS, k=[42]),
#         site_order = expand(rules.split_baypass_global_input_files.output.site_order, n=BAYPASS_SPLITS),
#     output:
#         c2_outliers = f"{ANALYSIS_DIR}/baypass/baypass_c2_outliers.txt",
#         c2_pval_hist = f"{ANALYSIS_DIR}/baypass/figures/baypass_c2_pval_hist.pdf",
#         c2_manhat_pdf = f"{ANALYSIS_DIR}/baypass/figures/baypass_c2_manhattan.pdf",
#         c2_manhat_png = f"{ANALYSIS_DIR}/baypass/figures/baypass_c2_manhattan.png"
#     params:
#         qval_cut = 0.05
#     conda: "../envs/baypass.yaml"
#     notebook:
#         "../notebooks/baypass_outlier_test.r.ipynb"

##############
#### POST ####
##############

rule baypass_done:
    input:
        # expand(rules.generate_windowed_c2_byCity.output, city=CITIES),
        # rules.fmd_and_omega_mat_pca.output,
        # rules.baypass_outlier_test.output,
        # expand(rules.baypass_coreModel_allSamples.output, city=CITIES, n=BAYPASS_SPLITS, k=[42]),
    output:
        f"{BAYPASS_DIR}/baypass.done"
    shell:
        """
        touch {output}
        """

