rule gridsearch_filter:
    input:
        cl=expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/absolute_PRE_down_sampling/clonality_results/{{project}}_{sample}_clonality.csv",sample=SAMPLES),
        rds=expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/absolute_PRE_down_sampling/relative_cn_rds/{{project}}_{sample}_{{bin}}kb_relSmoothedCN.rds",sample=SAMPLES)
    output:
        OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/{project}_fit_QC_predownsample.tsv"
    params:
        bin="{bin}",
        meta=config["samplesheet"],
        project="{project}",
        outdir=OUT_DIR,
        purity_cutoff=config["purity_cutoff"]
    threads: THREADS 
    script: 
        "../scripts/gridsearch_results_filtering_custom.R"
        
