#######################
#### IDENTIFY SNPS ####
#######################

rule angsd_snps_allSamples:
    """
    Identify SNPs across all samples using ANGSD
    """
    input:
        bams = rules.create_bam_list_allSamples_allSites.output,
        ref = rules.copy_ref.output,
        ref_idx = rules.samtools_index_ref.output
    output:
        gls = f'{ANGSD_DIR}/snps/allSamples/{{chrom}}/{{chrom}}_allSamples_snps.beagle.gz',
        mafs = f'{ANGSD_DIR}/snps/allSamples/{{chrom}}/{{chrom}}_allSamples_snps.mafs.gz',
        snp_stats = f'{ANGSD_DIR}/snps/allSamples/{{chrom}}/{{chrom}}_allSamples_snps.snpStat.gz',
        hwe = f'{ANGSD_DIR}/snps/allSamples/{{chrom}}/{{chrom}}_allSamples_snps.hwe.gz',
        pos = f'{ANGSD_DIR}/snps/allSamples/{{chrom}}/{{chrom}}_allSamples_snps.pos.gz',
        counts = f'{ANGSD_DIR}/snps/allSamples/{{chrom}}/{{chrom}}_allSamples_snps.counts.gz'
    log: f"{LOG_DIR}/angsd_snps_allSamples/{{chrom}}_angsd_snps.log"
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    params:
        max_depth = 5225, # Num samples x Mean coverage x 2
        min_ind = 1045, # 50% of Num samples
        maf = 0.05,
        out = f'{ANGSD_DIR}/snps/allSamples/{{chrom}}/{{chrom}}_allSamples_snps'
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 24000,
        runtime = 1440
    shell:
        """
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doMajorMinor 1 \
            -SNP_pval 1e-6 \
            -doMaf 1 \
            -minMaf {params.maf} \
            -doGlf 2 \
            -baq 2 \
            -ref {input.ref} \
            -doCounts 1 \
            -dumpCounts 4 \
            -setMinDepthInd 1 \
            -minInd {params.min_ind} \
            -setMaxDepth {params.max_depth} \
            -minQ 20 \
            -minMapQ 30 \
            -remove_bads 1 \
            -skipTriallelic 1 \
            -uniqueOnly 1 \
            -only_proper_pairs 1 \
            -dosnpstat 1 \
            -doHWE 1 \
            -r {wildcards.chrom} \
            -bam {input.bams} 2> {log}
        """

##########################################
#### EXCLUDE PUTATIVE PARALOGOUS SNPs ####
##########################################

rule create_ngsparalog_posfile:
    """
    Create position file of ngsParalog
    """
    input:
        rules.angsd_snps_allSamples.output.pos
    output:
        f"{PROGRAM_RESOURCE_DIR}/ngsparalog/{{chrom}}_posfile.txt"
    shell:
        """
        zcat {input} | tail -n +2 | cut -f1,2 > {output}
        """

rule ngsparalog:
    """
    Run ngsParalog in parallel across chromosomes
    """
    input:
        bams = rules.create_bam_list_allSamples_allSites.output,
        pos = rules.create_ngsparalog_posfile.output
    output:
        out = f"{NGSPARALOG_DIR}/{{chrom}}_ngsparalog.txt"
    log: f"{LOG_DIR}/ngsparalog/{{chrom}}_ngsparalog.log"
    container: "/home/santang3/singularity_containers/ngsparalog.sif"
    params:
        min_ind = 1045
    shell:
        """
        ( samtools mpileup -b {input.bams} \
            -l {input.pos} \
            -r {wildcards.chrom} \
            -q 0 -Q 0 --ff UNMAP,DUP |\
            ngsParalog calcLR \
                -infile - \
                -minQ 20 -minind {params.min_ind} -mincov 1 \
                -outfile {output} ) 2> {log}
        """

rule identify_paralogous_snps:
    """
    Identify SNPs in putatively paralogous alignments. Write sites files that exclude these
    """
    input:
        para = expand(rules.ngsparalog.output, chrom=CHROMOSOMES)
    output:
        manhat = f"{ANALYSIS_DIR}/ngsparalog/ngsparalog_manhattan.pdf",
        sites = expand(f"{PROGRAM_RESOURCE_DIR}/angsd_sites/{{chrom}}_filtered.sites", chrom=CHROMOSOMES)
    conda: "../envs/r.yaml"
    params:
        out = lambda w: f"{PROGRAM_RESOURCE_DIR}/angsd_sites"
    notebook:
        "../notebooks/identify_paralogous_snps.r.ipynb"

rule index_filtered_snps:
    """
    Index SNPs sites files for ANGSD
    """
    input:
        sites = lambda w: [x for x in rules.identify_paralogous_snps.output.sites if w.chrom in x]
    output:
        idx = f"{PROGRAM_RESOURCE_DIR}/angsd_sites/{{chrom}}_filtered.sites.idx",
        bin = f"{PROGRAM_RESOURCE_DIR}/angsd_sites/{{chrom}}_filtered.sites.bin"
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    shell:
        """
        angsd sites index {input}
        """

################################################
#### ANGSD ASSOCIATION BASED ON FREQUENCIES ####
################################################

rule create_angsd_asso_ybin_file:
    """
    Create ybin (i.e., phenotype file) for ANGSD association analysis
    """
    input:
        samples = config["samples"],
        bams = rules.create_bam_list_allSamples_allSites.output
    output:
        f"{PROGRAM_RESOURCE_DIR}/angsd_habitats.ybin"
    conda: "../envs/r.yaml"
    script:
        "../scripts/r/create_angsd_asso_ybin_file.R"

rule angsd_asso_freq:
    """
    Perform association analysis in ANGSD with binary variable (Urban vs. Rural)
    """
    input:
        ybin = rules.create_angsd_asso_ybin_file.output,
        bams = rules.create_bam_list_allSamples_allSites.output,
        ref = rules.copy_ref.output,
        ref_idx = rules.samtools_index_ref.output,
        sites = lambda w: [x for x in rules.identify_paralogous_snps.output.sites if w.chrom in x]
    output:
        asso = f'{ANGSD_DIR}/asso/allSamples/{{chrom}}/{{chrom}}_allSamples_freq.lrt0.gz'
    log: f"{LOG_DIR}/angsd_asso/{{chrom}}_angsd_asso_freq.log"
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    params:
        max_depth = 5225, # Num samples x Mean coverage x 2
        min_ind = 1045, # 50% of Num samples
        out = f'{ANGSD_DIR}/asso/allSamples/{{chrom}}/{{chrom}}_allSamples_freq'
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        runtime = 1440
    shell:
        """
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doMajorMinor 1 \
            -doMaf 1 \
            -baq 2 \
            -ref {input.ref} \
            -doCounts 1 \
            -setMinDepthInd 1 \
            -minInd {params.min_ind} \
            -setMaxDepth {params.max_depth} \
            -minQ 20 \
            -minMapQ 30 \
            -remove_bads 1 \
            -skipTriallelic 1 \
            -uniqueOnly 1 \
            -only_proper_pairs 1 \
            -doAsso 1 \
            -Pvalue 1 \
            -yBin {input.ybin} \
            -r {wildcards.chrom} \
            -sites {input.sites} \
            -bam {input.bams} 2> {log}
        """

####################################
#### ANGSD GENOTYPE LIKELIHOODS ####
####################################

rule angsd_gls_allSamples:
    """
    Estimate BEAGLE-formated genotype likelihoods for SNPs across all samples
    """
    input:
        bams = rules.create_bam_list_allSamples_allSites.output,
        ref = rules.copy_ref.output,
        ref_idx = rules.samtools_index_ref.output,
        sites = lambda w: [x for x in rules.identify_paralogous_snps.output.sites if w.chrom in x]
    output:
        gls = f'{ANGSD_DIR}/gls/allSamples/{{chrom}}/{{chrom}}_allSamples.beagle.gz',
    log: f"{LOG_DIR}/angsd_gl_allSamples/{{chrom}}_angsd_gl.log"
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    params:
        max_depth = 5225, # Num samples x Mean coverage x 2
        min_ind = 1045, # 50% of Num samples
        out = f'{ANGSD_DIR}/gls/allSamples/{{chrom}}/{{chrom}}_allSamples'
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        runtime = 1440
    shell:
        """
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doGlf 2 \
            -doMajorMinor 1 \
            -doCounts 1 \
            -setMinDepthInd 1 \
            -minInd {params.min_ind} \
            -setMaxDepth {params.max_depth} \
            -baq 2 \
            -ref {input.ref} \
            -minQ 20 \
            -minMapQ 30 \
            -remove_bads 1 \
            -skipTriallelic 1 \
            -uniqueOnly 1 \
            -only_proper_pairs 1 \
            -r {wildcards.chrom} \
            -sites {input.sites} \
            -bam {input.bams} 2> {log}
        """

##############################
#### POPULATION STRUCTURE ####
##############################

rule create_pos_file_for_ngsLD:
    input:
        rules.angsd_gls_allSamples.output.gls
    output:
        f'{PROGRAM_RESOURCE_DIR}/ngsld/{{chrom}}.pos'
    shell:
        """
        zcat {input} | tail -n +2 | cut -f1 | sed 's/_/\t/2' > {output}
        """
rule ngsld:
    input:
        pos = rules.create_pos_file_for_ngsLD.output,
        gls = rules.angsd_gls_allSamples.output.gls
    output:
        f'{NGSLD_DIR}/{{chrom}}.ld.gz'
    log: f'{LOG_DIR}/ngsld/{{chrom}}_ngsld.log'
    container: 'library://james-s-santangelo/ngsld/ngsld:1.2.0'
    threads: 4
    params:
        num_ind = 2090
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        runtime = 1440
    shell:
        """
        ( NUM_SITES=$(cat {input.pos} | wc -l) &&
          ngsLD --geno {input.gls} \
            --pos {input.pos} \
            --n_ind {params.num_ind} \
            --n_sites $NUM_SITES \
            --probs \
            --extend_out \
            --n_threads {threads} \
            --max_kb_dist 50 | gzip --best > {output} ) 2> {log}
        """

# rule angsd_asso_score:
#     """
#     Perform association analysis in ANGSD with binary variable (Urban vs. Rural)
#     """
#     input:
#         ybin = rules.create_angsd_asso_ybin_file.output,
#         bams = rules.create_bam_list_allSamples_allSites.output,
#         ref = rules.copy_ref.output,
#         ref_idx = rules.samtools_index_ref.output
#     output:
#         asso = f'{ANGSD_DIR}/asso/allSamples/{{chrom}}/{{chrom}}_allSamples_asso_score.lrt0.gz'
#     log: f"{LOG_DIR}/angsd_asso/{{chrom}}_angsd_asso_score.log"
#     container: 'library://james-s-santangelo/angsd/angsd:0.938'
#     params:
#         out = f'{ANGSD_DIR}/asso/allSamples/{{chrom}}/{{chrom}}_allSamples_asso_score'
#     threads: 8
#     resources:
#         mem_mb = lambda wildcards, attempt: attempt * 80000,
#         runtime = 1440
#     shell:
#         """
#         NUM_IND=$( wc -l < {input.bams} );
#         MIN_IND=$(( NUM_IND*50/100 ));
#         MAX_DEPTH=$(( NUM_IND*1*2 ));
#         angsd -GL 1 \
#             -out {params.out} \
#             -nThreads {threads} \
#             -doMajorMinor 1 \
#             -SNP_pval 1e-6 \
#             -doMaf 1 \
#             -minMaf 0.05 \
#             -baq 2 \
#             -ref {input.ref} \
#             -doCounts 1 \
#             -setMinDepthInd 1 \
#             -minInd $MIN_IND \
#             -setMaxDepth $MAX_DEPTH \
#             -minQ 20 \
#             -minMapQ 30 \
#             -doAsso 2 \
#             -Pvalue 1 \
#             -yBin {input.ybin} \
#             -doPost 1 \
#             -r {wildcards.chrom} \
#             -bam {input.bams} 2> {log}
#         """

# rule angsd_asso_lg:
#     """
#     Perform association analysis in ANGSD with binary variable (Urban vs. Rural)
#     """
#     input:
#         ybin = rules.create_angsd_asso_ybin_file.output,
#         bams = rules.create_bam_list_allSamples_allSites.output,
#         ref = rules.copy_ref.output,
#         ref_idx = rules.samtools_index_ref.output
#     output:
#         asso = f'{ANGSD_DIR}/asso/allSamples/{{chrom}}/{{chrom}}_allSamples_asso_lg.lrt0.gz'
#     log: f"{LOG_DIR}/angsd_asso/{{chrom}}_angsd_asso_lg.log"
#     container: 'library://james-s-santangelo/angsd/angsd:0.938'
#     params:
#         out = f'{ANGSD_DIR}/asso/allSamples/{{chrom}}/{{chrom}}_allSamples_asso_lg'
#     threads: 8
#     resources:
#         mem_mb = lambda wildcards, attempt: attempt * 80000,
#         runtime = 1440
#     shell:
#         """
#         NUM_IND=$( wc -l < {input.bams} );
#         MIN_IND=$(( NUM_IND*50/100 ));
#         MAX_DEPTH=$(( NUM_IND*1*2 ));
#         angsd -GL 1 \
#             -out {params.out} \
#             -nThreads {threads} \
#             -doMajorMinor 1 \
#             -SNP_pval 1e-6 \
#             -doMaf 1 \
#             -minMaf 0.05 \
#             -baq 2 \
#             -ref {input.ref} \
#             -doCounts 1 \
#             -setMinDepthInd 1 \
#             -minInd $MIN_IND \
#             -setMaxDepth $MAX_DEPTH \
#             -minQ 20 \
#             -minMapQ 30 \
#             -doAsso 4 \
#             -Pvalue 1 \
#             -yBin {input.ybin} \
#             -doPost 1 \
#             -r {wildcards.chrom} \
#             -bam {input.bams} 2> {log}
#         """

# rule analyze_angsd_asso:
#     input:
#         expand(rules.angsd_asso_freq.output, chrom=CHROMOSOMES[0]),
#         expand(rules.angsd_asso_lg.output, chrom=CHROMOSOMES[0]),
#         expand(rules.angsd_asso_score.output, chrom=CHROMOSOMES[0])
#     output:
#         "test.txt"
#     conda: "../envs/r.yaml"
#     notebook:
#         "../notebooks/analyze_angsd_asso.r.ipynb"

rule angsd_asso_done:
    input:
        expand(rules.angsd_asso_freq.output, chrom=CHROMOSOMES),
        expand(rules.ngsld.output, chrom=CHROMOSOMES)
    output:
        f"{ANGSD_DIR}/angsd_asso.done"
    shell:
        """
        touch {output}
        """
