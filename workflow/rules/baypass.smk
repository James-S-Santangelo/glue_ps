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
        perCity_geno = expand(f"{PROGRAM_RESOURCE_DIR}/baypass/allSites/{{city}}/{{city}}.geno", city=CITIES),
        perCity_cont = expand(f"{PROGRAM_RESOURCE_DIR}/baypass/allSites/{{city}}/{{city}}.cf", city=CITIES),
        site_order = f"{PROGRAM_RESOURCE_DIR}/baypass/allSites/site_order.txt",
        miss = f"{PROGRAM_RESOURCE_DIR}/baypass/missing.sites"
    conda: "../envs/baypass.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 60000,
        runtime = lambda wildcards, attempt: attempt * 120
    params:
        cities = CITIES,
        habitats = HABITATS,
        chrom = CHROMOSOMES
    notebook:
        "../notebooks/create_alleleCount_files_byCity.py.ipynb"


rule create_baypass_global_input_files:
    input:
        geno = expand(rules.create_alleleCount_files_byCity.output.perCity_geno, city=CITIES),
        cont = expand(rules.create_alleleCount_files_byCity.output.perCity_cont, city=CITIES)
    output:
        geno = f"{PROGRAM_RESOURCE_DIR}/baypass/allSites/allSamples.geno",
        cont = f"{PROGRAM_RESOURCE_DIR}/baypass/allSites/allSamples.cont",
    shell:
        """
        paste {input.geno} | tr '\t' ' ' > {output.geno}
        paste {input.cont} | tr '\t' ' ' > {output.cont}
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
    conda: "../envs/baypass.yaml"
    params:
        cities = CITIES,
        splits = BAYPASS_SPLITS
    notebook:
        "../notebooks/split_baypass_global_input_files.py.ipynb"

#################
#### BAYPASS ####
#################

def get_baypass_coreModel_input_geno(wildcards):
    all_geno_files = rules.split_baypass_global_input_files.output.as_geno
    geno = []
    for geno_file in all_geno_files:
        base = os.path.basename(geno_file)
        n = base.split(".")[0].split("_")[1]
        if n == wildcards.n:
            geno.append(geno_file)
    return geno

rule baypass_coreModel_allSamples:
    input:
        geno = get_baypass_coreModel_input_geno
    output:
        log = f"{BAYPASS_DIR}/coreModel_allSamples/{{n}}/allSamples_{{n}}_baypass.log",
        dic = f"{BAYPASS_DIR}/coreModel_allSamples/{{n}}/allSamples_{{n}}_DIC.out",
        omega_mat = f"{BAYPASS_DIR}/coreModel_allSamples/{{n}}/allSamples_{{n}}_mat_omega.out",
        beta_sum = f"{BAYPASS_DIR}/coreModel_allSamples/{{n}}/allSamples_{{n}}_summary_beta_params.out",
        omega_lda = f"{BAYPASS_DIR}/coreModel_allSamples/{{n}}/allSamples_{{n}}_summary_lda_omega.out",
        pif_sum = f"{BAYPASS_DIR}/coreModel_allSamples/{{n}}/allSamples_{{n}}_summary_pij.out",
        pi_xtx = f"{BAYPASS_DIR}/coreModel_allSamples/{{n}}/allSamples_{{n}}_summary_pi_xtx.out",
    log: f"{LOG_DIR}/baypass/coreModel_allSamples_{{n}}.log"
    container: "library://james-s-santangelo/baypass/baypass:2.41"
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        runtime = lambda wildcards, attempt: attempt * 720
    params:
        out_prefix = f"{BAYPASS_DIR}/coreModel_allSamples/allSamples_{{n}}"
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
        expand(rules.baypass_coreModel_allSamples.output, n=BAYPASS_SPLITS)
    output:
        f"{BAYPASS_DIR}/baypass.done"
    shell:
        """
        touch {output}
        """

