# Python functions used throughout snakemake workflow

def get_habitat_saf_files_allSites(wildcards):
    all_saf_files = expand(rules.angsd_saf_likelihood_byHabitat_allSites.output.saf_idx, habitat=HABITATS, chrom=wildcards.chrom, city=wildcards.city)
    first_hab = wildcards.hab_comb.split('_')[0]
    second_hab = wildcards.hab_comb.split('_')[1]
    saf1 = [x for x in all_saf_files if '{0}'.format(first_hab) in os.path.basename(x)]
    saf2 = [x for x in all_saf_files if '{0}'.format(second_hab) in os.path.basename(x)]
    return saf1 + saf2

def get_baypass_coreModel_input_geno(wildcards):
    all_geno_files = rules.split_baypass_global_input_files.output.as_geno
    geno = []
    for geno_file in all_geno_files:
        base = os.path.basename(geno_file)
        n = base.split(".")[0].split("_")[1]
        if n == wildcards.n:
            geno.append(geno_file)
    return geno

def get_baypass_coreModel_byCity_input(wildcards):
    all_geno_files = rules.split_baypass_global_input_files.output.perCity_geno
    all_cont_files = rules.create_alleleCount_files_byCity.output.perCity_cont
    all_pool_files = rules.create_alleleCount_files_byCity.output.perCity_pool
    geno = []
    for geno_file in all_geno_files:
        base = os.path.basename(geno_file)
        city = re.search("(.*)(?=_\\d+\\.geno)", base).group(1) 
        n = re.search("(\\d+)(?=\\.geno)", base).group(1)
        if n == wildcards.n and city == wildcards.city:
            geno.append(geno_file)
    cont = [f for f in all_cont_files if wildcards.city in os.path.basename(f)]
    pool = [f for f in all_pool_files if wildcards.city in os.path.basename(f)]
    return {"geno": geno, "cont": cont, "pool": pool}


