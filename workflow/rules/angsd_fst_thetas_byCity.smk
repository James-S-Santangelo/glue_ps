# Rules to estimate SAF/SFS by habitat using all sites. Used for detecting selective sweeps

###############
#### SETUP ####
###############

rule create_bam_list_allSamples_allSites:
    input:
        BAM_FILES_PATH
    output:
        f'{PROGRAM_RESOURCE_DIR}/bam_lists/allSamples_allSites_bams.list'
    log: LOG_DIR + '/create_bam_list/allSamples_allSites_bam_list.log'
    run:
        glue_bams = [bam for bam in glob.glob(f"{input}/final/*.bam") if not os.path.basename(bam).startswith('s_')]
        tor_bams = [bam for bam in glob.glob(f"{input}/toronto_bams/*.bam")]

        all_bams = glue_bams + tor_bams
        with open(output[0], 'w') as f:
            for bam in all_bams:
                sample = os.path.basename(bam).split('_merged')[0]
                if sample in SAMPLES:
                    f.write('{0}\n'.format(bam))


rule create_bam_list_byHabitat_allSites:
    """
    Create text file with paths to BAMS for urban, rural, and suburban samples. BAMs contain all reads genome-wide
    """
    input:
        bams = rules.create_bam_list_allSamples_allSites.output
    output:
        f'{PROGRAM_RESOURCE_DIR}/bam_lists/by_city/{{city}}/{{city}}_{{habitat}}_allSites_bams.list'
    params:
        samples = config['samples']
    run:
        df = pd.read_table(params.samples, sep = '\t')
        df_sub = df[(df['city'] == wildcards.city) & (df['site'] == wildcards.habitat)]
        samples_city_habitat = df_sub['sample'].tolist()
        bams = open(input[0], 'r').readlines()
        with open(output[0], 'w') as f:
            for bam in bams:
                sample = os.path.basename(bam).split('_merged')[0]
                if sample in samples_city_habitat:
                    f.write('{0}'.format(bam))

################################
#### SAF AND SFS ESTIMATION ####
################################

rule angsd_saf_likelihood_byHabitat_allSites:
    """
    Generate Site Allele Frequency (SAF) likelihood file for each habitat using ANGSD. Uses only 4fold sites.
    """
    input:
        bams = rules.create_bam_list_byHabitat_allSites.output,
        ref = rules.copy_ref.output,
        ref_idx = rules.samtools_index_ref.output
    output:
        saf = f'{ANGSD_DIR}/saf/by_city/{{city}}/{{habitat}}/{{chrom}}/{{chrom}}_{{habitat}}_allSites.saf.gz',
        saf_idx = f'{ANGSD_DIR}/saf/by_city/{{city}}/{{habitat}}/{{chrom}}/{{chrom}}_{{habitat}}_allSites.saf.idx',
        saf_pos = f'{ANGSD_DIR}/saf/by_city/{{city}}/{{habitat}}/{{chrom}}/{{chrom}}_{{habitat}}_allSites.saf.pos.gz'
    log: LOG_DIR + '/angsd_saf_likelihood_byHabitat_allSites/{city}/{city}_{habitat}_{chrom}_allSites_saf.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    params:
        out = f'{ANGSD_DIR}/saf/by_city/{{city}}/{{habitat}}/{{chrom}}/{{chrom}}_{{habitat}}_allSites'
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 10000,
        runtime = lambda wildcards, attempt: attempt * 180
    shell:
        """
        NUM_IND=$( wc -l < {input.bams} );
        MIN_IND=$(( NUM_IND*50/100 ));
        MAX_DEPTH=$(( NUM_IND*1*2 ));
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doMajorMinor 4 \
            -baq 2 \
            -ref {input.ref} \
            -doCounts 1 \
            -setMinDepthInd 1 \
            -minInd $MIN_IND \
            -setMaxDepth $MAX_DEPTH \
            -minQ 20 \
            -minMapQ 30 \
            -doSaf 1 \
            -anc {input.ref} \
            -r {wildcards.chrom} \
            -bam {input.bams} 2> {log}
        """
 
rule angsd_estimate_joint_habitat_sfs_allSites:
    """
    Estimated folded, two-dimensional urban-rural SFS for each city using realSFS. Uses all sites4
    """
    input:
        safs = get_habitat_saf_files_allSites
    output:
        f'{ANGSD_DIR}/sfs/{{city}}/2dsfs/allSites/{{chrom}}/{{city}}_{{chrom}}_allSites_{{hab_comb}}.2dsfs'
    log: LOG_DIR + '/angsd_estimate_habitat_2dsfs_allSites/{city}/{city}_{chrom}_allSites_{hab_comb}.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    threads: 12
    resources:
        mem_mb = lambda wildcards, attempt: 3000 if wildcards.city != 'Toronto' else 6000 ,
        runtime = lambda wildcards, attempt:  attempt * 360 if wildcards.city != 'Toronto' else attempt * 1440
    shell:
        """
        realSFS {input.safs} \
            -tole 1e-6 \
            -maxIter 30000 \
            -seed 42 \
            -fold 1 \
            -P {threads} > {output} 2> {log}
        """

rule angsd_estimate_sfs_byHabitat_allSites:
    """
    Estimate folded SFS separately for each habitat (i.e., 1D SFS) using realSFS. 
    """
    input:
        saf = rules.angsd_saf_likelihood_byHabitat_allSites.output.saf_idx
    output:
        f'{ANGSD_DIR}/sfs/{{city}}/1dsfs/allSites/{{chrom}}/{{city}}_{{chrom}}_allSites_{{habitat}}.sfs'
    log: LOG_DIR + '/angsd_estimate_sfs_byHabitat_allSites/{city}/{city}_{chrom}_allSites_{habitat}_sfs.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: 2000 if wildcards.city != 'Toronto' else 3000 ,
        runtime = lambda wildcards, attempt:  attempt * 180 if wildcards.city != 'Toronto' else attempt * 720
    shell:
        """
        realSFS {input.saf} \
            -P {threads} \
            -tole 1e-6 \
            -fold 1 \
            -maxIter 30000 \
            -seed 42 > {output} 2> {log}
        """

########################
### FST AND THETAS #####
########################

rule angsd_habitat_fst_index_allSites:
    """
    Estimate per-site alphas (numerator) and betas (denominator) for Fst estimation
    """
    input: 
        saf_idx = get_habitat_saf_files_allSites,
        joint_sfs = rules.angsd_estimate_joint_habitat_sfs_allSites.output
    output:
        fst = f'{ANGSD_DIR}/summary_stats/hudson_fst/allSites/{{city}}/{{chrom}}/{{city}}_{{chrom}}_allSites_{{hab_comb}}.fst.gz',
        idx = f'{ANGSD_DIR}/summary_stats/hudson_fst/allSites/{{city}}/{{chrom}}/{{city}}_{{chrom}}_allSites_{{hab_comb}}.fst.idx'
    log: LOG_DIR + '/angsd_habitat_fst_index_allSites/{city}_{chrom}_allSites_{hab_comb}_index.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    threads: 4
    resources:
        mem_mb = 4000,
        runtime = 120
    params:
        fstout = f'{ANGSD_DIR}/summary_stats/hudson_fst/allSites/{{city}}/{{chrom}}/{{city}}_{{chrom}}_allSites_{{hab_comb}}'
    shell:
        """
        realSFS fst index {input.saf_idx} \
            -sfs {input.joint_sfs} \
            -fold 1 \
            -P {threads} \
            -whichFst 1 \
            -fstout {params.fstout} 2> {log}
        """

rule angsd_fst_allSites_readable:
    """
    Create readable Fst files. Required due to format of realSFS fst index output files. 
    """
    input:
        rules.angsd_habitat_fst_index_allSites.output.idx
    output:
        f'{ANGSD_DIR}/summary_stats/hudson_fst/allSites/{{city}}/{{chrom}}/{{city}}_{{chrom}}_allSites_{{hab_comb}}_readable.fst'
    log: 'logs/angsd_fst_allSites_readable/{city}_{chrom}_{hab_comb}_readable_fst.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    shell:
        """
        realSFS fst print {input} > {output} 2> {log}
        """

rule angsd_estimate_thetas_byHabitat_allSites:
    """
    Generate per-site thetas in each habitat from 1DSFS
    """
    input:
        saf_idx = rules.angsd_saf_likelihood_byHabitat_allSites.output.saf_idx,
        sfs = rules.angsd_estimate_sfs_byHabitat_allSites.output
    output:
        idx = f'{ANGSD_DIR}/summary_stats/thetas/allSites/{{city}}/{{chrom}}/{{city}}_{{chrom}}_allSites_{{habitat}}.thetas.idx',
        thet = f'{ANGSD_DIR}/summary_stats/thetas/allSites/{{city}}/{{chrom}}/{{city}}_{{chrom}}_allSites_{{habitat}}.thetas.gz'
    log: LOG_DIR + '/angsd_estimate_thetas_byHabitat_allSites/{city}_{chrom}_allSites_{habitat}_thetas.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    threads: 4
    params:
        out = f'{ANGSD_DIR}/summary_stats/thetas/allSites/{{city}}/{{chrom}}/{{city}}_{{chrom}}_allSites_{{habitat}}'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        runtime = 60
    shell:
        """
        realSFS saf2theta {input.saf_idx} \
            -P {threads} \
            -fold 1 \
            -sfs {input.sfs} \
            -outname {params.out} 2> {log}
        """

rule angsd_thetas_allSites_readable:
    """
    Create readable Fst files. Required due to format of realSFS fst index output files. 
    """
    input:
        rules.angsd_estimate_thetas_byHabitat_allSites.output.idx
    output:
        f'{ANGSD_DIR}/summary_stats/thetas/allSites/{{city}}/{{chrom}}/{{city}}_{{chrom}}_allSites_{{habitat}}_readable.thetas'
    log: 'logs/angsd_thetas_allSites_readable/{city}_{chrom}_{habitat}_readable_thetas.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    shell:
        """
        thetaStat print {input} > {output} 2> {log}
        """

###########################
### WINDOWED ANALYSES #####
###########################

rule windowed_theta:
    input:
        rules.angsd_estimate_thetas_byHabitat_allSites.output.idx
    output:
        f"{ANGSD_DIR}/summary_stats/thetas/allSites/{{city}}/{{chrom}}/{{city}}_{{chrom}}_allSites_{{habitat}}_windowedThetas.gz.pestPG"
    log: LOG_DIR + '/windowed_theta/{city}_{chrom}_{habitat}_windowTheta.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    params:
        out = f"{ANGSD_DIR}/summary_stats/thetas/allSites/{{city}}/{{chrom}}/{{city}}_{{chrom}}_allSites_{{habitat}}_windowedThetas.gz",
        win = 20000,
        step = 20000
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        runtime = 60
    shell:
        """
        thetaStat do_stat {input} -win {params.win} -step {params.step} -outnames {params.out} 2> {log}
        """

rule windowed_fst:
    input:
        rules.angsd_habitat_fst_index_allSites.output.idx
    output:
        f"{ANGSD_DIR}/summary_stats/hudson_fst/allSites/{{city}}/{{chrom}}/{{city}}_{{chrom}}_allSites_{{hab_comb}}_windowed.fst"
    log: LOG_DIR + '/windowed_fst/{city}_{chrom}_{hab_comb}_windowedFst.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    params:
        win = 20000,
        step = 20000
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        runtime = 60
    shell:
        """
        realSFS fst stats2 {input} -win {params.win} -step {params.step} > {output} 2> {log}
        """

##############
#### POST ####
##############

rule angsd_fst_thetas_byCity_done:
    """
    Generate empty flag file signalling successful completion of SFS and summary stat for habitats
    """
    input:
        expand(rules.windowed_fst.output, city=CITIES, hab_comb=['urban_rural'], chrom=CHROMOSOMES), 
        expand(rules.windowed_theta.output, city=CITIES, habitat=HABITATS, chrom=CHROMOSOMES) 
    output:
        f'{ANGSD_DIR}/angsd_fsts_thetas_byCity.done'
    shell:
        """
        touch {output}
        """
