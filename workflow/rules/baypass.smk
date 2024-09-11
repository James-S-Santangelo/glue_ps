###############
#### SETUP ####
###############

rule create_alleleCount_files_byCity:
    input:
        perCity_mafs = lambda w: expand(rules.angsd_alleleCounts_byCity_byHabitat.output.mafs, city=CITIES, chrom=CHROMOSOMES, habitat=HABITATS),
        global_mafs = expand(rules.angsd_snp_af_allSamples.output.mafs, chrom=CHROMOSOMES),
        global_sites = expand(rules.extract_and_index_af_allSamples_sites.output.sites, chrom=CHROMOSOMES),
        pos = lambda w: expand(rules.angsd_alleleCounts_byCity_byHabitat.output.pos, city=CITIES, chrom=CHROMOSOMES, habitat=HABITATS),
        counts = lambda w: expand(rules.angsd_alleleCounts_byCity_byHabitat.output.counts, city=CITIES, chrom=CHROMOSOMES, habitat=HABITATS),
    output:
        perCity_geno = expand(f"{PROGRAM_RESOURCE_DIR}/baypass/{{city}}/{{city}}.geno", city=CITIES),
        perCity_cont = expand(f"{PROGRAM_RESOURCE_DIR}/baypass/{{city}}/{{city}}.cf", city=CITIES),
        site_order = f"{PROGRAM_RESOURCE_DIR}/baypass/site_order.txt",
        miss = f"{PROGRAM_RESOURCE_DIR}/baypass/missing.sites"
    conda: "../envs/baypass.yaml"
    params:
        cities = CITIES,
        habitats = HABITATS,
        chrom = CHROMOSOMES
    notebook:
        "../notebooks/create_baypass_input_files.py.ipynb"


rule create_baypass_global_input_files:
    input:
        geno = expand(rules.create_alleleCount_files_byCity.output.perCity_geno, city=CITIES),
        cont = expand(rules.create_alleleCount_files_byCity.output.perCity_cont, city=CITIES)
    output:
        geno = f"{PROGRAM_RESOURCE_DIR}/baypass/allSamples.geno",
        cont = f"{PROGRAM_RESOURCE_DIR}/baypass/allSamples.cont",
    shell:
        """
        paste {input.geno} | tr '\t' ' ' > {output.geno}
        paste {input.cont} | tr '\t' ' ' > {output.cont}
        """

#################
#### BAYPASS ####
#################

rule baypass_coreModel_allSamples:
    input:
        geno = rules.create_baypass_global_input_files.output.geno
    output:
        log = f"{BAYPASS_DIR}/coreModel_allSamples/allSamples_baypass.log",
        dic = f"{BAYPASS_DIR}/coreModel_allSamples/allSamples_DIC.out",
        omega_mat = f"{BAYPASS_DIR}/coreModel_allSamples/allSamples_mat_omega.out",
        beta_sum = f"{BAYPASS_DIR}/coreModel_allSamples/allSamples_summary_beta_params.out",
        omega_lda = f"{BAYPASS_DIR}/coreModel_allSamples/allSamples_summary_lda_omega.out",
        pif_sum = f"{BAYPASS_DIR}/coreModel_allSamples/allSamples_summary_pij.out",
        pi_xtx = f"{BAYPASS_DIR}/coreModel_allSamples/allSamples_summary_pi_xtx.out",
    log: f"{LOG_DIR}/baypass/coreModel_allSamples.log"
    container: "library://james-s-santangelo/baypass/baypass:2.41"
    threads: 12
    params:
        out_prefix = f"{BAYPASS_DIR}/coreModel_allSamples/allSamples"
    shell:
        """
        baypass -gfile {input.geno} \
            -outprefix {params.out_prefix} \
            -nthreads {threads} 2> {log}
        """


##############
#### POST ####
##############

rule baypass_done:
    input:
        rules.baypass_coreModel_allSamples.output 
    output:
        f"{BAYPASS_DIR}/baypass.done"
    shell:
        """
        touch {output}
        """

