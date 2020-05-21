# snakemake -s snake_VD.smk --jobs 5000 --latency-wait 30 --cluster-config cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources} -M {cluster.memory}  -o {cluster.output} -e {cluster.error}' --keep-going --rerun-incomplete
import numpy as np

configfile: "config_vd.yaml"

rule all:
    input:
        expand("{phenodir}/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Ysim_reg_{method}",
            phenodir=config["phenodir"],
            N=config["N"],
            seed=range(1,101),
            a=["{:.1f}".format(x) for x in np.append(0.05, np.arange(start=0.1,stop=1, step=0.1))],
            P=['50', '100', '1000'],
            h2=['0.3'],
            S=config["S"],
            model="modelnoiseFixedAndBggeneticFixedAndBg",
            method=['limmbo.csv', 'gemma.txt', 'sbat.txt']),
        expand("{phenodir}/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Ysim_reg_{method}",
            phenodir=config["phenodir"],
            N=config["N"],
            seed=range(101,201),
            a = "1",
            P=['10', '20', '30', '40', '50', '70', '100', '150', '200', '500', '1000'],
            h2=['0.3'],
            S=config["S"],
            model="modelnoiseFixedAndBggeneticFixedAndBg",
            method=['limmbo.csv', 'gemma.txt', 'sbat.txt']),
        expand("{phenodir}/Calibration/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Ysim_reg_{method}",
            phenodir=config["phenodir"],
            N='1000',
            seed=range(201,301),
            P=['10', '30', '100', '300'],
            h2=['0.3'],
            S="0",
            model="modelnoiseFixedAndBggeneticBgOnly",
            method=['limmbo.csv', 'gemma.txt', 'sbat.txt']),
        expand("{phenodir}/Scalability/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Ysim_reg_{method}",
            phenodir=config["phenodir"],
            N=['1000', '2000', '5000', '10000'],
            seed=range(201,301),
            P='100',
            h2=['0.3'],
            S="0",
            model="modelnoiseFixedAndBggeneticBgOnly",
            method=['limmbo.csv', 'gemma.txt', 'sbat.txt']),
        kinship=expand("{dir}/genotypes/relatedEU_nopopstructure/N{N}/{N}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_kinship_{method}",
            dir=config["simulatedir"],
            N=config["N"], A=config["paramGeno"]["nAnc"],
            B=config["paramGeno"]["nBlocks"],
            type=config["paramGeno"]["type"],
            seed=config["paramGeno"]["seedGeno"],
            method=['eigenvec_limmbo.csv', 'eigenvalue_limmbo.csv', 'eigenvec.txt', 'eigenvalue.txt'])

rule preparePhenotypes:
    input:
        pheno="{dir}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Ysim.csv",
        covs="{dir}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Covs.csv"
    output:
        "{dir}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Ysim_reg_gemma.txt",
        "{dir}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Ysim_reg_sbat.txt",
        "{dir}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Ysim_reg_limmbo.csv"
    shell:
        "Rscript 'preparePhenotypes.R' --pheno={input.pheno} --covs={input.covs}"


rule prepareKinship:
    input:
        kinship="{dir}/genotypes/relatedEU_nopopstructure/N{N}/{N}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_kinship_norm.csv"
    output:
        "{dir}/genotypes/relatedEU_nopopstructure/N{N}/{N}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_kinship_eigenvalue_limmbo.csv",
        "{dir}/genotypes/relatedEU_nopopstructure/N{N}/{N}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_kinship_eigenvec_limmbo.csv",
        "{dir}/genotypes/relatedEU_nopopstructure/N{N}/{N}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_kinship_eigenvalue.txt",
        "{dir}/genotypes/relatedEU_nopopstructure/N{N}/{N}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_kinship_eigenvec.txt"
    shell:
        "Rscript 'prepareKinship.R' --kinship={input.kinship}"

rule limmbo:
    input:
        kinship=expand("{dir}/genotypes/relatedEU_nopopstructure/N{N}/{N}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_kinship_norm.csv",
            dir=config["simulatedir"],
            N=config["N"], A=config["paramGeno"]["nAnc"],
            B=config["paramGeno"]["nBlocks"],
            type=config["paramGeno"]["type"],
            seed=config["paramGeno"]["seedGeno"]),
        pheno=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Ysim_reg_limmbo.csv",
            dir=config["simulatedir"])
    params:
        sp=config["params"]["sp"],
        min=config["params"]["min"],
        tr=config["params"]["tr"],
        cpus=config["params"]["cpu"]
    output:
        "{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Cg_fit_seed{seed}.csv",
        "{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Cn_fit_seed{seed}.csv"
    log:
        "{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/log/sim.log"
    #conda:
    #    "env/limmbo3.5.yaml"
    shell:
        "(runVarianceEstimation \
            -o {wildcards.dir}/vd/seedLiMMBo{wildcards.seed}/{wildcards.method} \
            -sp {params.sp} -seed {wildcards.seed} \
            -p {input.pheno} -k {input.kinship} \
            --minCooccurrence {params.min} \
            -cpus {params.cpus} --limmbo -v -t) 2> {log}"

rule limix:
    input:
        kinship=expand("{dir}/genotypes/relatedEU_nopopstructure/N{N}/{N}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_kinship_norm.csv",
            dir=config["simulatedir"],
            N=config["N"], A=config["paramGeno"]["nAnc"],
            B=config["paramGeno"]["nBlocks"],
            type=config["paramGeno"]["type"],
            seed=config["paramGeno"]["seedGeno"]),
        pheno=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Ysim_reg_limmbo.csv",
            dir=config["simulatedir"])
    params:
        tr=config["params"]["tr"],
    output:
        "{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Cg_REML_seed{seed}.csv",
        "{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Cn_REML_seed{seed}.csv",
    log:
        "{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/log/sim.log"
    #conda:
    #    "env/limmbo3.5.yaml"
    shell:
        "(runVarianceEstimation \
            -o {wildcards.dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed} \
            -seed {wildcards.seed} \
            -p {input.pheno} -k {input.kinship} \
            --reml -v -t) 2> {log}"
rule sbat:
    input:
        eigenvec=expand("{dir}/genotypes/relatedEU_nopopstructure/N{{N}}/{{N}}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_kinship_eigenvec.txt",
            dir=config["simulatedir"],
            A=config["paramGeno"]["nAnc"],
            B=config["paramGeno"]["nBlocks"],
            type=config["paramGeno"]["type"],
            seed=config["paramGeno"]["seedGeno"]),
        eigenvalue=expand("{dir}/genotypes/relatedEU_nopopstructure/N{N}/{N}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_kinship_eigenvalue.txt",
            dir=config["simulatedir"],
            N=config["N"], A=config["paramGeno"]["nAnc"],
            B=config["paramGeno"]["nBlocks"],
            type=config["paramGeno"]["type"],
            seed=config["paramGeno"]["seedGeno"]),
        pheno=expand("{dir}/phenotypes/Calibration/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Ysim_reg_sbat.txt",
            dir=config["simulatedir"])
    output:
        "{dir}/vd/Calibration/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/sbat.seed{seed}.C.mpmm.txt"
    log:
        "{dir}/vd/Calibration/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/log/sbat.seed{seed}.log"
    shell:
        "(sbat_fs --covariances-only \
              -p {input.pheno} \
              -ev {input.eigenvec} {input.eigenvalue} \
              -o {wildcards.dir}/vd/Calibration/Traits{wildcards.P}_samples{wildcards.N}_NrSNP{wildcards.S}_Cg{wildcards.h2}_{wildcards.model}/seed{wildcards.seed}/sbat.seed{wildcards.seed}) 2> {log}"

rule gemma:
    input:
        eigenvec=expand("{dir}/genotypes/relatedEU_nopopstructure/N{N}/{N}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_kinship_eigenvec.csv",
            dir=config["simulatedir"],
            N=config["N"], A=config["paramGeno"]["nAnc"],
            B=config["paramGeno"]["nBlocks"],
            type=config["paramGeno"]["type"],
            seed=config["paramGeno"]["seedGeno"]),
        eigenval=expand("{dir}/genotypes/relatedEU_nopopstructure/N{N}/{N}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_kinship_eigenvalue.csv",
            dir=config["simulatedir"],
            N=config["N"], A=config["paramGeno"]["nAnc"],
            B=config["paramGeno"]["nBlocks"],
            type=config["paramGeno"]["type"],
            seed=config["paramGeno"]["seedGeno"]),
        pheno=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Ysim_reg_gemma.txt",
            dir=config["simulatedir"])
    output:
        "{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/log/gemma.seed{seed}.log"
        "{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/sbat.seed{seed}.C.mpmm.txt",
        "{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/sbat.seed{seed}.D.mpmm.txt",
    params:
        num=lambda wildcards: print (' '.join(map(str, range(1,(wildcards.P+1)))))
    log:
    shell:
        """
        (gemma -p {input.pheno} \
           -d {input.eigenval} \
           -u {input.eigenvec} \
           -lmm 2 \
           -n {params.num} \
           -o {dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}) 2> {log}"
        

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
