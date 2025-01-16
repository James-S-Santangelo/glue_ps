# Python functions used throughout snakemake workflow

def get_habitat_saf_files_allSites(wildcards):
    all_saf_files = expand(rules.angsd_saf_likelihood_byHabitat_allSites.output.saf_idx, habitat=HABITATS, chrom=wildcards.chrom, city=wildcards.city)
    first_hab = wildcards.hab_comb.split('_')[0]
    second_hab = wildcards.hab_comb.split('_')[1]
    saf1 = [x for x in all_saf_files if '{0}'.format(first_hab) in os.path.basename(x)]
    saf2 = [x for x in all_saf_files if '{0}'.format(second_hab) in os.path.basename(x)]
    return saf1 + saf2

