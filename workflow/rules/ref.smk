# Rules for dealing with reference genome

rule copy_ref:
    input:
        REFERENCE_GENOME
    output:
        f"{REF_DIR}/TrR_v6_haploid_reference.fasta"
    conda: '../envs/ref.yaml'
    log: f"{LOG_DIR}/copy_ref/copy_ref.log"
    shell:
        """
        cp {input} {output}
        """

rule samtools_index_ref:
    input:
        rules.copy_ref.output
    output:
        f"{REF_DIR}/TrR_v6_haploid_reference.fasta.fai"
    conda: '../envs/ref.yaml'
    log: f"{LOG_DIR}/samtools_index_ref/samtools_index_ref.log"
    shell:
        """
        sleep 10
        samtools faidx {input} 2> {log}
        """
