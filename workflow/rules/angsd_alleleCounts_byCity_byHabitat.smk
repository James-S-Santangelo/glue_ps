# Rules to estimate allele counts for each habitat in each city.
# First estimates SNPs across global samples. Then counts allales at these SNPs in each habitat per city
# Will be used for running BayPass with binary contrast statistic

################################
#### GLOBAL SNP ESTIMATION ##### 
################################

rule angsd_snp_af_allSamples:
    """
    Generate allele frequencies at SNPs across all global samples
    """
    input:
        bams = rules.create_bam_list_allSamples_allSites.output,
        ref = rules.copy_ref.output,
        ref_idx = rules.samtools_index_ref.output
    output:
        mafs = f'{ANGSD_DIR}/maf/allSamples/{{chrom}}/{{chrom}}_allSamples_snps.mafs.gz',
    log: LOG_DIR + '/angsd_snp_af_allSamples/{chrom}_allSamples_snps.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    params:
        out = f'{ANGSD_DIR}/maf/allSamples/{{chrom}}/{{chrom}}_allSamples_snps'
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 10000,
        time = lambda wildcards, attempt: str(attempt * 3) + ":00:00" 
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
            -doMaf 1 \
            -SNP_pval 1e-6 \
            -setMinDepthInd 1 \
            -minInd $MIN_IND \
            -setMaxDepth $MAX_DEPTH \
            -minQ 20 \
            -minMapQ 30 \
            -anc {input.ref} \
            -r {wildcards.chrom} \
            -bam {input.bams} 2> {log}
        """

############################################
#### ALLELE COUNTS PER-CITY PER-HABITAT ####
############################################

rule extract_and_index_af_allSamples_sites:
    input:
        rules.angsd_snp_af_allSamples.output
    output:
        sites = f'{PROGRAM_RESOURCE_DIR}/angsd_sites/{{chrom}}/{{chrom}}_af_allSamples.sites',
        binary= f'{PROGRAM_RESOURCE_DIR}/angsd_sites/{{chrom}}/{{chrom}}_af_allSamples.sites.bin',
        idx = f'{PROGRAM_RESOURCE_DIR}/angsd_sites/{{chrom}}/{{chrom}}_af_allSamples.sites.idx'
    log: LOG_DIR + '/extract_and_index_af_allSamples_sites/{chrom}_angsd_index.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    shell:
        """
        ( zcat {input} | tail -n +2 | cut -f 1,2 > {output.sites} &&
        sleep 5
        angsd sites index {output.sites} ) 2> {log}
        """

rule angsd_alleleCounts_byCity_byHabitat:
    input:
        bams = rules.create_bam_list_byHabitat_allSites.output,
        sites = rules.extract_and_index_af_allSamples_sites.output.sites,
        ref = rules.copy_ref.output,
        ref_idx = rules.samtools_index_ref.output
    output:
        mafs = f'{ANGSD_DIR}/maf/byCity/{{city}}/{{chrom}}/{{city}}_{{habitat}}_{{chrom}}_snps.mafs.gz',
        pos = f'{ANGSD_DIR}/maf/byCity/{{city}}/{{chrom}}/{{city}}_{{habitat}}_{{chrom}}_snps.pos.gz',
        counts = f'{ANGSD_DIR}/maf/byCity/{{city}}/{{chrom}}/{{city}}_{{habitat}}_{{chrom}}_snps.counts.gz',
    log: LOG_DIR + '/angsd_alleleCounts_byCity_byHabitat/{city}/{chrom}_{habitat}_byCity_snps.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    params:
        out = f'{ANGSD_DIR}/maf/byCity/{{city}}/{{chrom}}/{{city}}_{{habitat}}_{{chrom}}_snps'
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 10000,
        time = lambda wildcards, attempt: str(attempt * 3) + ":00:00" 
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
            -dumpCounts 3 \
            -doMaf 1 \
            -minQ 20 \
            -minMapQ 30 \
            -sites {input.sites} \
            -anc {input.ref} \
            -r {wildcards.chrom} \
            -bam {input.bams} 2> {log}
        """
 
##############
#### POST ####
##############

rule angsd_alleleCounts_byCity_byHabitat_done:
    """
    Generate empty flag file signalling successful completion of allele counts by city 
    """
    input:
        expand(rules.angsd_alleleCounts_byCity_byHabitat.output, city=CITIES, habitat=HABITATS, chrom=CHROMOSOMES), 
    output:
        f'{ANGSD_DIR}/angsd_alleleCounts_byCity_byHabitat.done'
    shell:
        """
        touch {output}
        """
