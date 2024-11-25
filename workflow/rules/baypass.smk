###############
#### SETUP ####
###############

rule create_baypass_global_input_files:
    input:
        perCity_mafs = expand(rules.angsd_alleleCounts_byCity_byHabitat.output.mafs, city=CITIES, chrom=CHROMOSOMES, habitat=HABITATS),
        global_mafs = expand(rules.angsd_snp_af_allSamples.output.mafs, chrom=CHROMOSOMES),
        global_sites = expand(rules.extract_and_index_af_allSamples_sites.output.sites, chrom=CHROMOSOMES),
        pos = lambda w: expand(rules.angsd_alleleCounts_byCity_byHabitat.output.pos, city=CITIES, chrom=CHROMOSOMES, habitat=HABITATS),
        counts = lambda w: expand(rules.angsd_alleleCounts_byCity_byHabitat.output.counts, city=CITIES, chrom=CHROMOSOMES, habitat=HABITATS),
    output:
        as_geno = f"{PROGRAM_RESOURCE_DIR}/baypass/allSites/allSamples.geno",
        as_cont = f"{PROGRAM_RESOURCE_DIR}/baypass/allSites/allSamples.cf",
        site_order = f"{PROGRAM_RESOURCE_DIR}/baypass/allSites/site_order_allSamples.txt",
        miss = f"{PROGRAM_RESOURCE_DIR}/baypass/allSites/missing_allSamples.sites"
    conda: "../envs/baypass.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 40000,
        runtime = lambda wildcards, attempt: attempt * 60 
    params:
        maf_prefix = f"{ANGSD_DIR}/maf",
        cities = CITIES,
        habitats = HABITATS,
        chrom = CHROMOSOMES
    notebook:
        "../notebooks/create_alleleCount_files_byCity.py.ipynb"

rule split_baypass_global_input_files:
    input:
        as_geno = rules.create_baypass_global_input_files.output.as_geno,
        site_order = rules.create_baypass_global_input_files.output.site_order
    output:
        as_geno = expand(f"{PROGRAM_RESOURCE_DIR}/baypass/split_files/allSamples/splits/allSamples_{{n}}.geno", n=BAYPASS_SPLITS),
        site_order = expand(f"{PROGRAM_RESOURCE_DIR}/baypass/split_files/allSamples/site_order/site_order_{{n}}.txt", n=BAYPASS_SPLITS)
    conda: "../envs/baypass.yaml"
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

rule baypass_coreModel_allSamples:
    input:
        geno = get_baypass_coreModel_input_geno,
        cont = rules.create_baypass_global_input_files.output.as_cont,
    output:
        log = f"{BAYPASS_DIR}/coreModel_allSamples/seed{{k}}/split{{n}}/allSamples_seed{{k}}_split{{n}}_baypass.log",
        dic = f"{BAYPASS_DIR}/coreModel_allSamples/seed{{k}}/split{{n}}/allSamples_seed{{k}}_split{{n}}_DIC.out",
        omega_mat = f"{BAYPASS_DIR}/coreModel_allSamples/seed{{k}}/split{{n}}/allSamples_seed{{k}}_split{{n}}_mat_omega.out",
        beta_sum = f"{BAYPASS_DIR}/coreModel_allSamples/seed{{k}}/split{{n}}/allSamples_seed{{k}}_split{{n}}_summary_beta_params.out",
        omega_lda = f"{BAYPASS_DIR}/coreModel_allSamples/seed{{k}}/split{{n}}/allSamples_seed{{k}}_split{{n}}_summary_lda_omega.out",
        pij_sum = f"{BAYPASS_DIR}/coreModel_allSamples/seed{{k}}/split{{n}}/allSamples_seed{{k}}_split{{n}}_summary_pij.out",
        pi_xtx = f"{BAYPASS_DIR}/coreModel_allSamples/seed{{k}}/split{{n}}/allSamples_seed{{k}}_split{{n}}_summary_pi_xtx.out",
        cont_out = f"{BAYPASS_DIR}/coreModel_allSamples/seed{{k}}/split{{n}}/allSamples_seed{{k}}_split{{n}}_summary_contrast.out",
        bf_out = f"{BAYPASS_DIR}/coreModel_allSamples/seed{{k}}/split{{n}}/allSamples_seed{{k}}_split{{n}}_summary_betai_reg.out"
    log: f"{LOG_DIR}/baypass/coreModel_allSamples/seed{{k}}/seed{{k}}_split{{n}}.log"
    container: "library://james-s-santangelo/baypass/baypass:2.41"
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        runtime = 720 
    params:
        out_prefix = f"{BAYPASS_DIR}/coreModel_allSamples/seed{{k}}/split{{n}}/allSamples_seed{{k}}_split{{n}}"
    shell:
        """
        baypass -gfile {input.geno} \
            -outprefix {params.out_prefix} \
            -seed {wildcards.k} \
            -contrastfile {input.cont} \
            -efile {input.cont} \
            -nthreads {threads} 2> {log}
        """

# rule baypass_coreModel_byCity:
#     input:
#         unpack(get_baypass_coreModel_byCity_input)
#     output:
#         log = f"{BAYPASS_DIR}/coreModel_byCity/{{city}}/split{{n}}/{{city}}_split{{n}}_baypass.log",
#         dic = f"{BAYPASS_DIR}/coreModel_byCity/{{city}}/split{{n}}/{{city}}_split{{n}}_DIC.out",
#         omega_mat = f"{BAYPASS_DIR}/coreModel_byCity/{{city}}/split{{n}}/{{city}}_split{{n}}_mat_omega.out",
#         beta_sum = f"{BAYPASS_DIR}/coreModel_byCity/{{city}}/split{{n}}/{{city}}_split{{n}}_summary_beta_params.out",
#         omega_lda = f"{BAYPASS_DIR}/coreModel_byCity/{{city}}/split{{n}}/{{city}}_split{{n}}_summary_lda_omega.out",
#         pif_sum = f"{BAYPASS_DIR}/coreModel_byCity/{{city}}/split{{n}}/{{city}}_split{{n}}_summary_pij.out",
#         pi_xtx = f"{BAYPASS_DIR}/coreModel_byCity/{{city}}/split{{n}}/{{city}}_split{{n}}_summary_pi_xtx.out",
#         cont_out = f"{BAYPASS_DIR}/coreModel_byCity/{{city}}/split{{n}}/{{city}}_split{{n}}_summary_contrast.out",
#         bf_out = f"{BAYPASS_DIR}/coreModel_byCity/{{city}}/split{{n}}/{{city}}_split{{n}}_summary_betai_reg.out"
#     log: f"{LOG_DIR}/baypass/coreModel_byCity/{{city}}/{{city}}_split{{n}}.log"
#     container: "library://james-s-santangelo/baypass/baypass:2.41"
#     threads: 4
#     resources:
#         mem_mb = lambda wildcards, attempt: attempt * 2000,
#         runtime = lambda wildcards, attempt: attempt * 240
#     params:
#         out_prefix = f"{BAYPASS_DIR}/coreModel_byCity/{{city}}/split{{n}}/{{city}}_split{{n}}"
#     shell:
#         """
#         baypass -gfile {input.geno} \
#             -outprefix {params.out_prefix} \
#             -seed 1 \
#             -contrastfile {input.cont} \
#             -efile {input.cont} \
#             -nthreads {threads} 2> {log}
#         """

##################
#### ANALYSES ####
##################

rule simulate_baypass:
    input:
        omega = expand(rules.baypass_coreModel_allSamples.output.omega_mat, n=37, k=3700),
        betas = expand(rules.baypass_coreModel_allSamples.output.beta_sum, n=37, k=3700),
        geno = lambda w: [x for x in rules.split_baypass_global_input_files.output.as_geno if "37" in x]
    output:
        simu_out = f"{PROGRAM_RESOURCE_DIR}/baypass/simulated.geno"
    conda: "../envs/baypass.yaml"
    notebook:
        "../notebooks/simulate_baypass.r.ipynb"
        
rule baypass_coreModel_simulated:
    input:
        geno = rules.simulate_baypass.output.simu_out,
        cont = rules.create_baypass_global_input_files.output.as_cont,
    output:
        log = f"{BAYPASS_DIR}/coreModel_simulated/simulated_baypass.log",
        dic = f"{BAYPASS_DIR}/coreModel_simulated/simulated_DIC.out",
        omega_mat = f"{BAYPASS_DIR}/coreModel_simulated/simulated_mat_omega.out",
        beta_sum = f"{BAYPASS_DIR}/coreModel_simulated/simulated_summary_beta_params.out",
        omega_lda = f"{BAYPASS_DIR}/coreModel_simulated/simulated_summary_lda_omega.out",
        pij_sum = f"{BAYPASS_DIR}/coreModel_simulated/simulated_summary_pij.out",
        pi_xtx = f"{BAYPASS_DIR}/coreModel_simulated/simulated_summary_pi_xtx.out",
        cont_out = f"{BAYPASS_DIR}/coreModel_simulated/simulated_summary_contrast.out",
        bf_out = f"{BAYPASS_DIR}/coreModel_simulated/simulated_summary_betai_reg.out"
    log: f"{LOG_DIR}/baypass/coreModel_simulated/coreModel_simulated.log"
    container: "library://james-s-santangelo/baypass/baypass:2.41"
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        runtime = 720 
    params:
        out_prefix = f"{BAYPASS_DIR}/coreModel_simulated/simulated"
    shell:
        """
        baypass -gfile {input.geno} \
            -outprefix {params.out_prefix} \
            -seed 42 \
            -contrastfile {input.cont} \
            -efile {input.cont} \
            -nthreads {threads} 2> {log}
        """

rule fmd_and_omega_mat_pca:
    input:
        obs_omega_mat = expand(rules.baypass_coreModel_allSamples.output.omega_mat, n=BAYPASS_SPLITS, k=[42, 420, 3700]),
        obs_beta_sum = expand(rules.baypass_coreModel_allSamples.output.beta_sum, n=BAYPASS_SPLITS, k=[42, 420, 3700]),
        sim_omega_mat = expand(rules.baypass_coreModel_simulated.output.omega_mat),
        sim_beta_sum = expand(rules.baypass_coreModel_simulated.output.beta_sum)
    output:
        fmd_box = f"{ANALYSIS_DIR}/baypass/figures/fmd_boxplot.pdf",
        fmd_sum = f"{ANALYSIS_DIR}/baypass/tables/fmd_summary.txt",
        alpha_box = f"{ANALYSIS_DIR}/baypass/figures/alpha_boxplot.pdf",
        beta_box = f"{ANALYSIS_DIR}/baypass/figures/beta_boxplot.pdf",
        obs_sim_stats = f"{ANALYSIS_DIR}/baypass/tables/observe_vs_simulated_statistics.txt",
        obs_rand_pca = f"{ANALYSIS_DIR}/baypass/figures/pca_by_continent_and_habitat_observed.pdf",
        sim_pca = f"{ANALYSIS_DIR}/baypass/figures/pca_by_continent_and_habitat_simulated.pdf",
    conda: "../envs/baypass.yaml"
    params:
        cities = CITIES,
        habitats = HABITATS,
    notebook:
        "../notebooks/fmd_and_omega_mat_pca.r.ipynb"

# rule generate_windowed_c2_byCity:
#     input:
#         c2 = lambda w: expand(rules.baypass_coreModel_byCity.output.cont_out, city=w.city, n=BAYPASS_SPLITS),
#         site_order = expand(rules.split_baypass_global_input_files.output.site_order, n=BAYPASS_SPLITS)
#     output:
#         win_c2 = f"{ANALYSIS_DIR}/baypass/windowed_c2/{{city}}_windowed_c2.txt"
#     conda: "../envs/baypass.yaml"
#     script:
#         "../scripts/r/generate_windowed_c2_byCity.R"

# rule baypass_outlier_test:
#     input:
#         cont_out = expand(rules.baypass_coreModel_allSamples.output.cont_out, n=BAYPASS_SPLITS, k=[1]),
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
        expand(rules.baypass_coreModel_allSamples.output, city=CITIES, n=BAYPASS_SPLITS, k=[42, 420, 3700]),
        rules.fmd_and_omega_mat_pca.output
        # expand(rules.baypass_coreModel_byCity.output, city=CITIES, n=BAYPASS_SPLITS)
    output:
        f"{BAYPASS_DIR}/baypass.done"
    shell:
        """
        touch {output}
        """

