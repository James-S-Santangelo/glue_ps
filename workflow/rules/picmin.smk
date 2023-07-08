rule install_picmin_dependencies:
    output:
       f"{ANALYSIS_DIR}/picmin/dependencies.done"
    conda: '../envs/picmin.yaml'
    shell:
        """
        Rscript -e "install.packages('poolr', repos = 'http://cran.us.r-project.org')"
        Rscript -e "devtools::install_github('TBooker/PicMin', force = T)"
        touch {output}
        """

rule generate_picmin_data:
    input:
        flag = rules.install_picmin_dependencies.output,
        fst = expand(rules.windowed_fst.output, city=CITIES, hab_comb=['urban_rural'], chrom=CHROMOSOMES), 
    output:
        fst = f"{ANALYSIS_DIR}/picmin/csv/all_windowed_fst.csv",
        hist = f"{ANALYSIS_DIR}/picmin/figures/missing_data_histogram.pdf"
    conda: '../envs/picmin.yaml'
    script:
        "../scripts/r/generate_picmin_data.R"

rule run_picmin:
    input:
        flag = rules.install_picmin_dependencies.output,
        fst = rules.generate_picmin_data.output.fst
    output:
        f"{ANALYSIS_DIR}/picmin/csv/outliers/outliers_{{n}}_cities.csv",
        f"{ANALYSIS_DIR}/picmin/figures/manhattan/picmin_manhattan_plot_{{n}}_cities.pdf",
        f"{ANALYSIS_DIR}/picmin/figures/outlier_histograms/outliers_histogram_{{n}}_cities.pdf"
    conda: '../envs/picmin.yaml'
    script:
        "../scripts/r/run_picmin.R"

rule picmin_done:
    input:
        expand(rules.run_picmin.output, n=[i for i in range(10, 27)])
    output:
        f"{ANALYSIS_DIR}/picmin/picmin.done"
    shell:
        """
        touch {output}
        """
