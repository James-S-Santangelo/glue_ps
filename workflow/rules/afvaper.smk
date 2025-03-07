rule generate_afvaper_input_files:
    input:
        afs = lambda w: expand(
            rules.angsd_alleleCounts_byCity_byHabitat.output.mafs, 
            city = CITIES, 
            habitat = HABITATS,
            chrom = w.chrom)
    output:
        f"{AFVAPER_DIR}/af_matrices/{{chrom}}.txt"
    conda: "../envs/r.yaml"
    notebook:
        "../notebooks/generate_afvaper_input_files.r.ipynb"

rule install_afvaper:
    output:
        f"{AFVAPER_DIR}/afvaper_install.done"
    conda: "../envs/r.yaml"
    shell:
        """
        R -e 'remotes::install_github("JimWhiting91/afvaper")' &&
        R -e 'library(afvaper)' &&
        touch {output}
        """

rule run_afvaper:
    input:
        install = rules.install_afvaper.output,
        fai = rules.samtools_index_ref.output,
        afs = rules.generate_afvaper_input_files.output
    output:
        f"{AFVAPER_DIR}/rds_files/{{chrom}}_afvaper.rds"
    conda: "../envs/r.yaml"
    script:
        "../scripts/r/run_afvaper.R"

rule afvaper_done:
    input:
        expand(rules.run_afvaper.output, chrom=CHROMOSOMES)
    output:
        f"{AFVAPER_DIR}/afvaper.done"
    shell:
        """
        touch {output}
        """
