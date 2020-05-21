# snakemake -s scalability.smk --jobs 5000 --latency-wait 30 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources} -M {cluster.memory}  -o {cluster.output} -e {cluster.error}' --keep-going --rerun-incomplete
import numpy as np

configfile: "config/config_scalability.yaml"

rule all:
    input:
        expand("{dir}/vd/Scalability/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/{type}",
            dir=config["dir"],
            N=['1000', '2000', '5000', '10000'],
            seed='201',
            P='100',
            h2=['0.3'],
            S="0",
            model="modelnoiseFixedAndBggeneticBgOnly",
            type=['sbat.C.mpmm.txt', 'Cg_gemma.csv'])

rule scalability_limmbo:
    input:
        eigenvec=expand("{dir}/phenotypes/Scalability/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/kinship_eigenvec_limmbo.csv",
            dir=config["simulatedir"]),
        eigenvalue=expand("{dir}/phenotypes/Scalability/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/kinship_eigenvalue_limmbo.csv",
            dir=config["simulatedir"]),
        pheno=expand("{dir}/phenotypes/Scalability/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Ysim_reg_limmbo.csv",
            dir=config["simulatedir"])
    params:
        sp=lambda wildcards: '5' if wildcards.P == 10 else '10',
        min=config["params"]["min"],
        tr=config["params"]["tr"],
        cpus=config["params"]["cpu"]
    output:
        "{dir}/vd/Scalability/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Cg_fit.csv",
        "{dir}/vd/Scalability/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Cn_fit.csv"
    conda:
        "envs/limmbo1.0.0-dev.yaml"
    shell:
        "(runVarianceEstimation \
            -o {wildcards.dir}/vd/Scalability/Traits{wildcards.P}_samples{wildcards.N}_NrSNP{wildcards.S}_Cg{wildcards.h2}_{wildcards.model}/seed{wildcards.seed} \
            -sp {params.sp} -seed {wildcards.seed} \
            -p {input.pheno} \
            -eval_k {input.eigenvalue} \
            -evec_k {input.eigenvec} \
            --minCooccurrence {params.min} \
            -cpus {params.cpus} --limmbo -v -t) 2> {log}"

rule scalability_limix:
    input:
        eigenvec=expand("{dir}/phenotypes/Scalability/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/kinship_eigenvec_limmbo.csv",
            dir=config["simulatedir"]),
        eigenvalue=expand("{dir}/phenotypes/Scalability/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/kinship_eigenvalue_limmbo.csv",
            dir=config["simulatedir"]),
        pheno=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Ysim_reg_limmbo.csv",
            dir=config["simulatedir"])
    params:
        tr=config["params"]["tr"],
    output:
        "{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Cg_REML_seed{seed}.csv",
        "{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Cn_REML_seed{seed}.csv",
    log:
        "{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/log/sim.log"
    conda:
        "envs/limmbo1.0.0-dev.yaml"
    shell:
        "(runVarianceEstimation \
            -o {wildcards.dir}/vd/Scalability/Traits{wildcards.P}_samples{wildcards.N}_NrSNP{wildcards.S}_Cg{wildcards.h2}_{wildcards.model}/seed{wildcards.seed} \
            -seed {wildcards.seed} \
            -p {input.pheno} \
            -eval_k {input.eigenvalue} \
            -evec_k {input.eigenvec} \
            --reml -v -t) 2> {log}"

rule scalability_sbat:
    input:
        eigenvec=expand("{dir}/phenotypes/Scalability/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/kinship_eigenvec.txt",
            dir=config["simulatedir"]),
        eigenvalue=expand("{dir}/phenotypes/Scalability/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/kinship_eigenvalue.txt",
            dir=config["simulatedir"]),
        pheno=expand("{dir}/phenotypes/Scalability/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Ysim_reg_sbat.txt",
            dir=config["simulatedir"])
    output:
        C="{dir}/vd/Scalability/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/sbat.C.mpmm.txt",
        D="{dir}/vd/Scalability/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/sbat.D.mpmm.txt",
        sbatlog="{dir}/vd/Scalability/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/sbat.log",
        time="{dir}/vd/Scalability/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/VD_time_sbat.csv"
    shell:
        """
        (sbat_fs --covariances-only \
            -p {input.pheno} \
            -ev {input.eigenvec} {input.eigenvalue} \
            --threads 1 \
            -o {wildcards.dir}/vd/Scalability/Traits{wildcards.P}_samples{wildcards.N}_NrSNP{wildcards.S}_Cg{wildcards.h2}_{wildcards.model}/seed{wildcards.seed}/sbat) 1>{output.sbatlog}
        grep 'Running time:' {output.sbatlog} | cut -d ":" -f 2 | cut -d " " -f 2 > {output.time}
        """

rule scalability_gemma:
    input:
        eigenvec=expand("{dir}/phenotypes/Scalability/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/kinship_eigenvec.txt",
            dir=config["simulatedir"]),
        eigenvalue=expand("{dir}/phenotypes/Scalability/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/kinship_eigenvalue.txt",
            dir=config["simulatedir"]),
        pheno=expand("{dir}/phenotypes/Scalability/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Ysim_reg_gemma.txt",
            dir=config["simulatedir"]),
        geno=expand("{dir}/genotypes/relatedEU_nopopstructure/N{{N}}/{{N}}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_chr1_maf{maf}_gemma.csv",
            dir=config["simulatedir"],
            A=config["paramGeno"]["nAnc"],
            B=config["paramGeno"]["nBlocks"],
            type=config["paramGeno"]["type"],
            maf=config["paramGeno"]["maf"],
            seed=config["paramGeno"]["seedGeno"])
    output:
        gemmalog="{dir}/vd/Scalability/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/gemma.log",
        time="{dir}/vd/Scalability/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/VD_time_gemma.csv",
        ve="{dir}/vd/Scalability/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Cn_gemma.csv",
        vg="{dir}/vd/Scalability/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Cg_gemma.csv"
    params:
        num=lambda wildcards: ' '.join(map(str, range(1, (int(wildcards.P) + 1))))
    log:
    shell:
        """
        (gemma -p {input.pheno} \
           -g {input.geno} \
           -d {input.eigenvalue} \
           -u {input.eigenvec} \
           -lmm 2 \
           -n {params.num} \
           -o {wildcards.dir}/vd/Scalability/Traits{wildcards.P}_samples{wildcards.N}_NrSNP{wildcards.S}_Cg{wildcards.h2}_{wildcards.model}/seed{wildcards.seed}) 2> {output.gemmalog}
        GemmaVD-parser.py -f {output.log} -o {wildcards.dir}/vd/Scalability/Traits{wildcards.P}_samples{wildcards.N}_NrSNP{wildcards.S}_Cg{wildcards.h2}_{wildcards.model}/seed{wildcards.seed}
        """

rule mvlmm:
    input:
        kinship=expand("{dir}/genotypes/relatedEU_nopopstructure/N{N}/{N}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_kinship_norm.csv",
            dir=config["simulatedir"],
            N=config["N"], A=config["paramGeno"]["nAnc"],
            B=config["paramGeno"]["nBlocks"],
            type=config["paramGeno"]["type"],
            seed=config["paramGeno"]["seedGeno"]),
        pheno=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Ysim.csv",
            dir=config["simulatedir"])
    output:
        "{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/sbat.seed{seed}.C.mpmm.txt",
        "{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/sbat.seed{seed}.D.mpmm.txt",
    log:
        "{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/log/gemma.log"
    shell:
        "(python /homes/hannah/software/anaconda2/envs/limmbo-paper-py2/bin/pylmmGWAS.py  -p [filename] \
               -k [filename] \
               -n [num] \
               -o [prefix]) 2> {log}"
