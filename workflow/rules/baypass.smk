###############
#### SETUP ####
###############

rule create_baypass_input_files_byCity:
    input:
        perCity_mafs = lambda w: expand(rules.angsd_alleleCounts_byCity_byHabitat.output.mafs, city=w.city, chrom=CHROMOSOMES, habitat=HABITATS),
        global_mafs = expand(rules.angsd_snp_af_allSamples.output.mafs, chrom=CHROMOSOMES),
        global_sites = expand(rules.extract_and_index_af_allSamples_sites.output.sites, chrom=CHROMOSOMES),
        pos = lambda w: expand(rules.angsd_alleleCounts_byCity_byHabitat.output.pos, city=w.city, chrom=CHROMOSOMES, habitat=HABITATS),
        counts = lambda w: expand(rules.angsd_alleleCounts_byCity_byHabitat.output.counts, city=w.city, chrom=CHROMOSOMES, habitat=HABITATS),
    output:
        perCity_geno = f"{PROGRAM_RESOURCE_DIR}/baypass/{{city}}/{{city}}.geno",
        perCity_cont = f"{PROGRAM_RESOURCE_DIR}/baypass/{{city}}/{{city}}.cf"
    conda: "../envs/baypass.yaml"
    params:
        cities = CITIES,
        habitats = HABITATS,
        chrom = CHROMOSOMES
    notebook:
        "notebooks/create_baypass_input_files.py.ipynb"

rule create_baypass_global_input_files:
    input:
        geno = expand(rules.create_baypass_input_files_byCity.output.perCity_geno, city=CITIES),
        cont = expand(rules.create_baypass_input_files_byCity.output.perCity_cont, city=CITIES)
    output:
        geno = f"{PROGRAM_RESOURCE_DIR}/baypass/allSamples.geno",
        cont = f"{PROGRAM_RESOURCE_DIR}/baypass/allSamples.cf"
    shell:
        """
        paste {input.geno} | tr '\t' ' ' > {output.geno}
        paste {input.cont} | tr '\t' ' ' > {output.cont}
        """

#################
#### BAYPASS ####
#################



##############
#### POST ####
##############

rule baypass_done:
    input:
        rules.create_baypass_global_input_files.output 
    output:
        f"{BAYPASS_DIR}/baypass.done"
    shell:
        """
        touch {output}
        """

