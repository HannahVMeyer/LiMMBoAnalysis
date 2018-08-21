# snakemake -s prepare.smk --jobs 5000 --latency-wait 30 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources} -M {cluster.memory}  -o {cluster.output} -e {cluster.error}' --keep-going --rerun-incomplete
import numpy as np

configfile: "config/config_prepare.yaml"

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
            N=['1000', '10000'], A=config["paramGeno"]["nAnc"],
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
        "Rscript 'prepare/preparePhenotypes.R' --pheno={input.pheno} --covs={input.covs}"


rule prepareKinship:
    input:
        kinship="{dir}/genotypes/relatedEU_nopopstructure/N{N}/{N}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_kinship_norm.csv"
    output:
        "{dir}/genotypes/relatedEU_nopopstructure/N{N}/{N}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_kinship_eigenvalue_limmbo.csv",
        "{dir}/genotypes/relatedEU_nopopstructure/N{N}/{N}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_kinship_eigenvec_limmbo.csv",
        "{dir}/genotypes/relatedEU_nopopstructure/N{N}/{N}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_kinship_eigenvalue.txt",
        "{dir}/genotypes/relatedEU_nopopstructure/N{N}/{N}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_kinship_eigenvec.txt"
    shell:
        "Rscript 'prepare/prepareKinship.R' --kinship={input.kinship}"
