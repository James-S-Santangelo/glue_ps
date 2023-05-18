# Python functions used throughout snakemake workflow

def run_fast_scandir(dir, sample):
    """
    Helper function to perform fast recursive search
    """
    subfolders, files = [], []

    for f in os.scandir(dir):
        if f.is_dir():
            subfolders.append(f.path)
        if f.is_file():
            if re.findall(r'^{0}_'.format(sample), f.name):
                files.append(f.path)
    for dir in list(subfolders):
        sf, f = run_fast_scandir(dir, sample)
        subfolders.extend(sf)
        files.extend(f) 
    return subfolders, sorted(files)

def get_raw_reads(wildcards):
    """
    Recursively search for forward and reverse reads for sample
    """
    folders, bams = run_fast_scandir(BAM_FILES, wildcards.sample)
    bam = reads[0]
    return bam

def get_files_for_saf_estimation_byHabitat(wildcards):
    ref = rules.unzip_reference.output
    bams = expand(rules.create_bam_list_byHabitat.output, habitat = wildcards.habitat, site=wildcards.site)
    sites = rules.convert_sites_for_angsd.output
    idx = rules.angsd_index_degenerate_sites.output
    chroms = config['chromosomes']
    return { 'bams' : bams, 'ref' : ref, 'sites' : sites, 'idx' : idx, 'chroms' : chroms }

def get_habitat_saf_files(wildcards):
    all_saf_files = expand(rules.angsd_saf_likelihood_byHabitat.output.saf_idx, habitat=HABITATS, site=wildcards.site)
    first_hab = wildcards.hab_comb.split('_')[0]
    second_hab = wildcards.hab_comb.split('_')[1]
    saf1 = [x for x in all_saf_files if '{0}'.format(first_hab) in os.path.basename(x)]
    saf2 = [x for x in all_saf_files if '{0}'.format(second_hab) in os.path.basename(x)]
    return saf1 + saf2

def get_habitat_saf_files_allSites(wildcards):
    all_saf_files = expand(rules.angsd_saf_likelihood_byHabitat_allSites.output.saf_idx, habitat=HABITATS, chrom=wildcards.chrom, city=wildcards.city)
    first_hab = wildcards.hab_comb.split('_')[0]
    second_hab = wildcards.hab_comb.split('_')[1]
    saf1 = [x for x in all_saf_files if '{0}'.format(first_hab) in os.path.basename(x)]
    saf2 = [x for x in all_saf_files if '{0}'.format(second_hab) in os.path.basename(x)]
    return saf1 + saf2

def get_files_for_saf_estimation_byPopulation(wildcards):
    ref = rules.unzip_reference.output
    bams = expand(rules.create_bam_list_byPop_multiInd.output, popu=wildcards.popu, site=wildcards.site)
    sites = rules.convert_sites_for_angsd.output
    idx = rules.angsd_index_degenerate_sites.output
    chroms = config['chromosomes']
    return { 'bams' : bams, 'ref' : ref, 'sites' : sites, 'idx' : idx, 'chroms' : chroms }

def get_population_saf_files(wildcards):
    all_saf_files = expand(rules.angsd_saf_likelihood_byPopulation.output.saf_idx, popu=POPS_MULTI_IND, site='4fold')
    pop1 = wildcards.pop_comb.split('_')[0]
    pop2 = wildcards.pop_comb.split('_')[1]
    saf1 = [x for x in all_saf_files if os.path.basename(x).startswith('{0}_'.format(pop1))]
    saf2 = [x for x in all_saf_files if os.path.basename(x).startswith('{0}_'.format(pop2))]
    return saf1 + saf2 

def get_whatshap_phase_input(wildcards):
    ref = REFERENCE_GENOME
    vcf = rules.bcftools_split_samples.output
    bam = expand(rules.samtools_markdup.output.bam, sample = wildcards.sample)
    return { 'ref' : ref, 'vcf' : vcf, 'bam' : bam }

def get_dadi_sfs_input_files(wildcards):
    hab1 = wildcards.hab_comb.split('_')[0]
    hab2 = wildcards.hab_comb.split('_')[1]
    saf_files = expand(rules.angsd_saf_likelihood_byHabitat.output.saf_idx, habitat=HABITATS, site='4fold')  
    sfs_files = expand(rules.angsd_estimate_sfs_byHabitat.output, habitat=HABITATS, site='4fold') 
    saf_urban = [x for x in saf_files if '{0}'.format(hab1) in os.path.basename(x)]
    saf_rural = [x for x in saf_files if '{0}'.format(hab2) in os.path.basename(x)]
    sfs_urban = [x for x in sfs_files if '{0}'.format(hab1) in os.path.basename(x)]
    sfs_rural = [x for x in sfs_files if '{0}'.format(hab2) in os.path.basename(x)]
    ref = rules.unzip_reference.output
    return { 'saf_urban' : saf_urban , 'saf_rural' : saf_rural, 'sfs_urban' : sfs_urban, 'sfs_rural' : sfs_rural, 'ref' : ref }

def selscan_xpnsl_input(wildcards):
    if 'Urban' in wildcards.hab_comb:
        vcf = expand(rules.bcftools_splitVCF_byHabitat.output.vcf, chrom=wildcards.chrom, habitat='Urban')
    elif 'Suburban' in wildcards.hab_comb:
        vcf = expand(rules.bcftools_splitVCF_byHabitat.output.vcf, chrom=wildcards.chrom, habitat='Suburban')
    vcf_ref = expand(rules.bcftools_splitVCF_byHabitat.output.vcf, chrom=wildcards.chrom, habitat='Rural')
    return { 'vcf' : vcf, 'vcf_ref' : vcf_ref } 

def xpclr_input(wildcards):
    if 'Urban' in wildcards.hab_comb:
        pop1s = expand(rules.samples_byHabitat.output, chrom=wildcards.chrom, habitat='Urban')
    elif 'Suburban' in wildcards.hab_comb:
        pop1s = expand(rules.samples_byHabitat.output, chrom=wildcards.chrom, habitat='Suburban')
    pop2s = expand(rules.samples_byHabitat.output, chrom=wildcards.chrom, habitat='Rural')
    vcf = expand(rules.shapeit_phase.output.vcf, chrom=wildcards.chrom)
    return { 'vcf' : vcf, 'pop1s' : pop1s, 'pop2s' : pop2s } 


