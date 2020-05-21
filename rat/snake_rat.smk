#snakemake -s snake_rat.smk --use-conda --jobs 5000 --latency-wait 30 --cluster-config cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources} -M {cluster.memory}  -o {cluster.output} -e {cluster.error}' --keep-going --rerun-incomplete ~/data/LiMMBo/rat/association/phenix/lmm_mt_betavalue_chr{1..21}.csv


configfile: "config_rat.yaml"

rule all:
    input:
        expand("{dir}/association/{method}/lmm_mt_betavalue_chr{chr}.csv",
            dir=config['directory'],
            method=['phenix', 'mvn', 'mice'],
            chr=range(1,21)
            ),
        expand("{dir}/association/{method}/lmm_mt_pvalue_chr{chr}.csv",
            dir=config['directory'],
            method=['phenix', 'mvn', 'mice'],
            chr=range(1,21)
            )

rule genotypes:
    input:
        hf5file="{dir}/rawdata/arrayexpress/HS.hdf5"
    params:
        dir="{dir}/rawdata/arrayexpress",
        name="HS_rats",
        estimateKin="True"
    output:
        "{dir}/rawdata/arrayexpress/HS_rats_kinship.csv"
    log:
        "{dir}/rawdata/arrayexpress/log/genotypes.log"
    script:
        "preprocessing/genotypes/01_ratdata-format-h5f.py"

rule limmbo:
    input:
        pheno="{dir}/processeddata/phenotypes_{method}_reg.csv",
        kinship="{dir}/rawdata/arrayexpress/HS_rats_kinship_norm.csv"
    params:
        sp=lambda wildcards: config["params"]["sp"],
        min=lambda wildcards: config["params"]["min"],
        tr=lambda wildcards: config["params"]["tr"],
        cpus=lambda wildcards: config["params"]["cpu"]
    output:
        "{dir}/vd/seedLiMMBo{seed}/{method}/Cg_fit_seed{seed}.csv",
        "{dir}/vd/seedLiMMBo{seed}/{method}/Cn_fit_seed{seed}.csv"
    #conda:
    #    "env/limmbo3.5.yaml"
    log:
        "{dir}/vd/seedLiMMBo{seed}/{method}/log/limmbo.log"
    shell:
        "(runVarianceEstimation \
        -o {wildcards.dir}/vd/seedLiMMBo{wildcards.seed}/{wildcards.method} \
            -sp {params.sp} -seed {wildcards.seed} \
            -p {input.pheno} -k {input.kinship} --minCooccurrence {params.min} \
            -tr {params.tr} -cpus {params.cpus} --limmbo -v -t) 2> {log}"

rule multitrait_limmbo:
    input:
        pheno="{dir}/processeddata/phenotypes_{method}_reg.csv",
        geno="{dir}/rawdata/arrayexpress/HS_rats_chr{chr}_maf0.05.csv",
        kinship="{dir}/rawdata/arrayexpress/HS_rats_kinship_norm.csv",
        cg=expand("{{dir}}/vd/seedLiMMBo{seed}/{{method}}/Cg_fit_seed{seed}.csv",
            seed= config["seed"]),
        cn=expand("{{dir}}/vd/seedLiMMBo{seed}/{{method}}/Cn_fit_seed{seed}.csv",
            seed= config["seed"])
    params:
        tr=lambda wildcards: config["params"]["tr"],
    output:
        "{dir}/association/{method}/lmm_mt_betavalue_chr{chr}.csv",
        "{dir}/association/{method}/lmm_mt_pvalue_chr{chr}.csv",
    #conda:
    #    "env/limmbo3.5.yaml"
    log:
        "{dir}/association/{method}/log/lmm_mt_chr{chr}.log"
    shell:
        "(runAssociation -mt -lmm \
            -p {input.pheno} -k {input.kinship} \
            -g {input.geno} \
            -o {wildcards.dir}/association/{wildcards.method} \
            -n chr{wildcards.chr} \
            -cg {input.cg} -cn {input.cn}\
            -tr {params.tr} --plot -v ) 2> {log}"

rule singletrait:
    input:
        pheno="{dir}/processeddata/phenotypes_{method}_reg.csv",
        geno="{dir}/rawdata/arrayexpress/HS_rats_chr{chr}_maf0.05.csv",
        kinship="{dir}/rawdata/arrayexpress/HS_rats_kinship_norm.csv",
    params:
        tr=lambda wildcards: config["params"]["tr"],
    output:
        "{dir}/association/{method}/lmm_st_betavalue_chr{chr}.csv",
        "{dir}/association/{method}/lmm_st_pvalue_chr{chr}.csv",
    #conda:
    #    "env/limmbo3.5.yaml"
    log:
        "{dir}/association/{method}/log/lmm_st_chr{chr}.log"
    shell:
        "(runAssociation -st -lmm \
            -p {input.pheno} -k {input.kinship} \
            -g {input.geno} \
            -o {wildcards.dir}/association/{wildcards.method} \
            -n chr{wildcards.chr} \
            -tr {params.tr} --plot -v ) 2> {log}"
