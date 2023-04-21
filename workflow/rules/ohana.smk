# Rules to run OHANA for detection of selection

###############
#### SETUP ####
###############

rule concat_habitat_bamLists_withinCities:
    """
    Concatenate urban and rural sample BAM lists within cities. Generates a single file with
    the paths to all of the BAM files for samples within a city
    """
    input:
        lambda w: expand(rules.create_bam_list_byHabitat_allSites.output, city=w.city, habitat=HABITATS)
    output:
        f'{PROGRAM_RESOURCE_DIR}/bam_lists/by_city/{{city}}/{{city}}_bams.list'
    shell:
        """
        cat {input} > {output} 
        """

######################
#### ESTIMATE GLs ####
######################

rule gls_allSamples:
    """
    Estimate BEAGLE-formated genotype likelihoods for SNPs across all samples
    """
    input:
        bams = rules.create_bam_list_allSamples_allSites.output,
        ref = rules.copy_ref.output,
        ref_idx = rules.samtools_index_ref.output
    output:
        gls = temp(f'{ANGSD_DIR}/gls/allSamples/{{chrom}}/{{chrom}}_allSamples.beagle.gz'),
        mafs = temp(f'{ANGSD_DIR}/gls/allSamples/{{chrom}}/{{chrom}}_allSamples.mafs.gz')
    log: f"{LOG_DIR}/angsd_gl_allSamples/{{chrom}}_angsd_gl.log"
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    params:
        out = f'{ANGSD_DIR}/gls/allSamples/{{chrom}}/{{chrom}}_allSamples'
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 80000,
        time = '24:00:00'
    shell:
        """
        NUM_IND=$( wc -l < {input.bams} );
        MIN_IND=$(( NUM_IND*50/100 ));
        MAX_DEPTH=$(( NUM_IND*1*2 ));   
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doGlf 2 \
            -doMajorMinor 1 \
            -SNP_pval 1e-6 \
            -doMaf 1 \
            -doCounts 1 \
            -setMinDepthInd 1 \
            -setMaxDepth $MAX_DEPTH \
            -baq 2 \
            -ref {input.ref} \
            -minInd $MIN_IND \
            -minQ 20 \
            -minMapQ 30 \
            -r {wildcards.chrom} \
            -bam {input.bams} 2> {log}
        """

rule concat_angsd_gl_allSamples:
    """
    Concatenated GLs from all 16 chromosomes into single file. Done separately for each site type.
    """
    input:
        expand(rules.gls_allSamples.output.gls, chrom=CHROMOSOMES)
    output:
        f'{ANGSD_DIR}/gls/allSamples/allChroms_allSamples.beagle.gz'
    log: f'{LOG_DIR}/concat_angsd_gl/allSamples_concat.log'
    conda: '../envs/ref.yaml'
    resources:
        time = '03:00:00',
        mem_mb = 4000
    shell:
        """
        first=1
        for f in {input}; do
            if [ "$first"  ]; then
                zcat "$f"
                first=
            else
                zcat "$f"| tail -n +2
            fi
        done | bgzip -c > {output} 2> {log}
        """

rule concat_angsd_mafs_allSamples:
    """
    Concatenate MAF files for each of 16 chromosomes into single file. Done separately for each site type.
    """
    input:
        expand(rules.gls_allSamples.output.mafs, chrom=CHROMOSOMES)
    output:
        f'{ANGSD_DIR}/gls/allSamples/allChroms_allSamples.mafs.gz'
    log: f'{LOG_DIR}/concat_angsd_mafs/allSamples_concat.log'
    conda: '../envs/ref.yaml'
    shell:
        """
        first=1
        for f in {input}; do
            if [ "$first"  ]; then
                zcat "$f"
                first=
            else
                zcat "$f"| tail -n +2
            fi
        done | bgzip -c > {output} 2> {log}
        """

rule gls_byCity:
    """
    Estimate BEAGLE-formatted genotype likelihoods for all samples within a particular city
    """
    input:
        bams = rules.concat_habitat_bamLists_withinCities.output,
        ref = rules.copy_ref.output
    output:
        gls = f'{ANGSD_DIR}/gls/by_city/{{city}}/{{city}}.beagle.gz',
        mafs = f'{ANGSD_DIR}/gls/by_city/{{city}}/{{city}}.mafs.gz'
    log: f'{LOG_DIR}/angsd_gl_byCity_beagle/{{city}}_beagleGL.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    params:
        out = f'{ANGSD_DIR}/gls/by_city/{{city}}/{{city}}',
        chroms = config['chromosomes']
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 12000,
        time = lambda wildcards, attempt: str(attempt * 24) + ":00:00" 
    shell:
        """
        NUM_IND=$( wc -l < {input.bams} );
        MIN_IND=$(( NUM_IND * 50/100 ));
        MAX_DEPTH=$(( NUM_IND*1*2 ));   
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doGlf 2 \
            -doMajorMinor 1 \
            -SNP_pval 1e-6 \
            -doMaf 1 \
            -doCounts 1 \
            -setMinDepthInd 1 \
            -baq 2 \
            -ref {input.ref} \
            -minInd $MIN_IND \
            -setMaxDepth $MAX_DEPTH \
            -minQ 20 \
            -minMapQ 30 \
            -rf {params.chroms} \
            -bam {input.bams} 2> {log}
        """

###############
#### OHANA ####
###############

rule convert_bglTolgm:
    input:
        rules.gls_allSamples.output.gls
    output:
        f"{OHANA_DIR}/bgl2lgm/{{chrom}}_allSamples.lgm"
    log: f"{LOG_DIR}/bgl2lgm/bgl2lgm_{{chrom}}_alLSamples.log"
    container: 'library://james-s-santangelo/ohana/ohana:latest'
    shell:
        """
        zcat {input} | convert bgl2lgm > {output} 2> {log}
        """

rule sample_lgm_sites:
    input:
        rules.convert_bglTolgm.output
    output:
         f"{OHANA_DIR}/lgm_downsampled/{{chrom}}_allSamples.lgm"
    log: f"{LOG_DIR}/sample_lgm_sites/{{chrom}}_sample-sites.log"
    container: 'library://james-s-santangelo/ohana/ohana:latest'
    params:
        percent = 5
    shell:
        """
        sample-sites.py {input} {params.percent} {output} 2> {log} 
        """

##############
#### POST ####
##############

rule ohana_done:
    input:
        expand(rules.sample_lgm_sites.output, chrom=CHROMOSOMES),
        rules.concat_angsd_gl_allSamples.output,
        rules.concat_angsd_mafs_allSamples.output,
        expand(rules.gls_byCity.output, city=CITIES)
    output:
        f'{OHANA_DIR}/ohana.done'
    shell:
        """
        touch {output}
        """

