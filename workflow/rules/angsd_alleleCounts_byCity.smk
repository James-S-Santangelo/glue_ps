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
 
##############
#### POST ####
##############

rule angsd_alleleCounts_byCity_done:
    """
    Generate empty flag file signalling successful completion of allele counts by city 
    """
    input:
        expand(rules.angsd_snp_af_allSamples.output, chrom=CHROMOSOMES), 
    output:
        f'{ANGSD_DIR}/angsd_alleleCounts_byCity.done'
    shell:
        """
        touch {output}
        """
