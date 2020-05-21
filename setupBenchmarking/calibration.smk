# snakemake -s calibration.smk --jobs 5000 --latency-wait 30 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources} -M {cluster.memory}  -o {cluster.output} -e {cluster.error}' --keep-going --rerun-incomplete
import numpy as np

configfile: "config/config_calibration.yaml"

rule all:
    input:
        expand("{dir}/vd/Calibration/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/{type}",
            dir=config["dir"],
            N=['5000'],
            seed=range(301,401),
            P=['10', '30'],
            h2=['0.3'],
            S="0",
            model="modelnoiseFixedAndBggeneticBgOnly",
            type=['Cg_REML.csv']),
            #type=['Cg_REML.csv', 'Cg_fit.csv', 'sbat.C.mpmm.txt', 'Cg_gemma.csv'])
        expand("{dir}/vd/Calibration/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/{type}",
            dir=config["dir"],
            N=['5000'],
            seed=range(301,401),
            P=['10', '30', '100', '300', '1000'],
            h2=['0.3'],
            S="0",
            model="modelnoiseFixedAndBggeneticBgOnly",
            type=['Cg_fit.csv', 'sbat.C.mpmm.txt'])


rule calibration_limmbo:
    input:
        eigenvec=expand("{dir}/phenotypes/Calibration/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Kinship_eigenvec_limmbo.csv",
            dir=config["simulatedir"]),
        eigenvalue=expand("{dir}/phenotypes/Calibration/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Kinship_eigenval_limmbo.csv",
            dir=config["simulatedir"]),
        pheno=expand("{dir}/phenotypes/Calibration/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Ysim_reg_limmbo.csv",
            dir=config["simulatedir"])
    params:
        sp=lambda wildcards: '5' if wildcards.P == 10 else '10',
        min=config["params"]["min"],
        tr=config["params"]["tr"],
        cpus=config["params"]["cpu"]
    output:
        "{dir}/vd/Calibration/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Cg_fit.csv",
        "{dir}/vd/Calibration/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Cn_fit.csv"
    log:
        "{dir}/vd/Calibration/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/log/vd.log"
    conda:
        "envs/limmbo1.0.0-dev.yaml"
    shell:
        "(runVarianceEstimation \
            -o {wildcards.dir}/vd/Calibration/Traits{wildcards.P}_samples{wildcards.N}_NrSNP{wildcards.S}_Cg{wildcards.h2}_{wildcards.model}/seed{wildcards.seed} \
            -sp {params.sp} -seed {wildcards.seed} \
            -p {input.pheno} \
            -eval_k {input.eigenvalue} \
            -evec_k {input.eigenvec} \
            --minCooccurrence {params.min} \
            -cpus {params.cpus} --limmbo -v -t) 2> {log}"

rule calibration_limix:
    input:
        eigenvec=expand("{dir}/phenotypes/Calibration/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Kinship_eigenvec_limmbo.csv",
            dir=config["simulatedir"]),
        eigenvalue=expand("{dir}/phenotypes/Calibration/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Kinship_eigenval_limmbo.csv",
            dir=config["simulatedir"]),
        pheno=expand("{dir}/phenotypes/Calibration/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Ysim_reg_limmbo.csv",
            dir=config["simulatedir"])
    params:
        tr=config["params"]["tr"],
    output:
        "{dir}/vd/Calibration/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Cg_REML.csv",
        "{dir}/vd/Calibration/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Cn_REML.csv",
    log:
        "{dir}/vd/Calibration/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/log/vd.log"
    conda:
         "envs/limmbo1.0.0-dev.yaml"
    shell:
        "(runVarianceEstimation \
            -o {wildcards.dir}/vd/Calibration/Traits{wildcards.P}_samples{wildcards.N}_NrSNP{wildcards.S}_Cg{wildcards.h2}_{wildcards.model}/seed{wildcards.seed} \
            -seed {wildcards.seed} \
            -eval_k {input.eigenvalue} \
            -evec_k {input.eigenvec} \
            -p {input.pheno} \
            --reml -v -t) 2> {log}"

rule calibration_sbat:
    input:
        eigenvec=expand("{dir}/phenotypes/Calibration/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Kinship_eigenvec_gemma.txt",
            dir=config["simulatedir"]),
        eigenvalue=expand("{dir}/phenotypes/Calibration/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Kinship_eigenval_gemma.txt",
            dir=config["simulatedir"]),
        pheno=expand("{dir}/phenotypes/Calibration/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Ysim_reg_sbat.txt",
            dir=config["simulatedir"])
    output:
        C="{dir}/vd/Calibration/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/sbat.C.mpmm.txt",
        D="{dir}/vd/Calibration/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/sbat.D.mpmm.txt",
        sbatlog="{dir}/vd/Calibration/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/sbat.log",
        time="{dir}/vd/Calibration/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/VD_time_sbat.csv"
    shell:
        """
        (sbat_fs --covariances-only \
            -p {input.pheno} \
            -ev {input.eigenvec} {input.eigenvalue} \
            -o {wildcards.dir}/vd/Calibration/Traits{wildcards.P}_samples{wildcards.N}_NrSNP{wildcards.S}_Cg{wildcards.h2}_{wildcards.model}/seed{wildcards.seed}/sbat) 1>{output.sbatlog}
        grep 'Running time:' {output.sbatlog} | cut -d ":" -f 2 | cut -d " " -f 2 > {output.time}
        """

rule calibration_gemma:
    input:
        eigenvec=expand("{dir}/phenotypes/Calibration/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Kinship_eigenvec_gemma.txt",
            dir=config["simulatedir"]),
        eigenval=expand("{dir}/phenotypes/Calibration/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Kinship_eigenval_gemma.txt",
            dir=config["simulatedir"]),
        pheno=expand("{dir}/phenotypes/Calibration/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Ysim_reg_gemma.txt",
            dir=config["simulatedir"]),
        geno=expand("{dir}/genotypes/relatedEU_nopopstructure/N{{N}}/{{N}}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_chr1_maf{maf}_gemma.csv",
            dir=config["simulatedir"],
            A=config["paramGeno"]["nAnc"],
            B=config["paramGeno"]["nBlocks"],
            type=config["paramGeno"]["type"],
            maf=config["paramGeno"]["maf"],
            seed=config["paramGeno"]["seedGeno"])
    output:
        gemmalog="{dir}/vd/Calibration/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/gemma.log",
        time="{dir}/vd/Calibration/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/VD_time_gemma.csv",
        ve="{dir}/vd/Calibration/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Cn_gemma.csv",
        vg="{dir}/vd/Calibration/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Cg_gemma.csv"
    params:
        num=lambda wildcards: ' '.join(map(str, range(1, (int(wildcards.P)+1))))
    shell:
        """
        (gemma0.96 -p {input.pheno} \
           -g <(head -n 1 {input.geno}) \
           -d {input.eigenval} \
           -u {input.eigenvec} \
           -lmm 2 \
           -n {params.num} \
           -o {wildcards.dir}/vd/Calibration/Traits{wildcards.P}_samples{wildcards.N}_NrSNP{wildcards.S}_Cg{wildcards.h2}_{wildcards.model}/seed{wildcards.seed}) 2> {output.gemmalog}
        CalibrationAndScalability/GemmaVD-parser.py -f {output.gemmalog} \
           -o {wildcards.dir}/vd/Calibration/Traits{wildcards.P}_samples{wildcards.N}_NrSNP{wildcards.S}_Cg{wildcards.h2}_{wildcards.model}/seed{wildcards.seed}
        """
