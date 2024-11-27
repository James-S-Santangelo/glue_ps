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
        thetas = expand(rules.windowed_theta.output, city=CITIES, habitat=['urban', 'rural'], chrom=CHROMOSOMES), 
    output:
        stats = f"{ANALYSIS_DIR}/picmin/csv/all_windowed_stats.csv",
        hist_fst = f"{ANALYSIS_DIR}/picmin/figures/missing_data_histogram_fst.pdf",
        hist_tp = f"{ANALYSIS_DIR}/picmin/figures/missing_data_histogram_tp.pdf",
        hist_td = f"{ANALYSIS_DIR}/picmin/figures/missing_data_histogram_td.pdf"
    conda: '../envs/picmin.yaml'
    script:
        "../scripts/r/generate_picmin_data.R"

rule run_picmin:
    input:
        flag = rules.install_picmin_dependencies.output,
        stats = rules.generate_picmin_data.output.stats
    output:
        picmin = f"{ANALYSIS_DIR}/picmin/csv/picmin/{{stat}}/picmin_{{stat}}_{{n}}_cities.csv",
        manhat = f"{ANALYSIS_DIR}/picmin/figures/manhattan/{{stat}}/picmin_manhattan_plot_{{stat}}_{{n}}_cities.pdf",
        hist = f"{ANALYSIS_DIR}/picmin/figures/outlier_histograms/{{stat}}/outliers_histogram_{{stat}}_{{n}}_cities.pdf"
    conda: '../envs/picmin.yaml'
    params:
        qval_cut = 0.05
    script:
        "../scripts/r/run_picmin.R"

rule pool_picmin:
    input:
        results = lambda w: expand(rules.run_picmin.output.picmin, stat=w.stat, n=[i for i in range(10, 27)])
    output:
        picmin = f"{ANALYSIS_DIR}/picmin/csv/picmin/{{stat}}/picmin_{{stat}}_pooled.csv",
        manhat = f"{ANALYSIS_DIR}/picmin/figures/manhattan/{{stat}}/picmin_manhattan_plot_{{stat}}_pooled.pdf",
        hist = f"{ANALYSIS_DIR}/picmin/figures/outlier_histograms/{{stat}}/outliers_histogram_{{stat}}_pooled.pdf"
    conda: '../envs/picmin.yaml'
    params:
        qval_cut = 0.05
    script:
        "../scripts/r/pool_picmin.R"

rule picmin_done:
    input:
        expand(rules.pool_picmin.output, stat=["fst", "td", "tp"]),
    output:
        f"{ANALYSIS_DIR}/picmin/picmin.done"
    shell:
        """
        touch {output}
        """
