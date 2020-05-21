# snakemake -s prepare.smk --jobs 5000 --latency-wait 30 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources} -M {cluster.memory}  -o {cluster.output} -e {cluster.error}' --keep-going --rerun-incomplete
import numpy as np

configfile: "config/config_prepare.yaml"

rule all:
    input:
        expand("{phenodir}/Calibration/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Ysim_reg_{method}",
            phenodir=config["phenodir"],
            N=['1000', '5000'],
            seed=range(301,401),
            P=['10', '30', '100', '300', '1000'],
            h2=['0.3'],
            S="0",
            model="modelnoiseFixedAndBggeneticBgOnly",
            method=['limmbo.csv', 'gemma.txt', 'sbat.txt']),

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
        kinship="{dir}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Kinship.csv"
    output:
        "{dir}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Kinship_eigenvalue_limmbo.csv",
        "{dir}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Kinship_eigenvec_limmbo.csv",
        "{dir}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Kinship_eigenvalue.txt",
        "{dir}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Kinship_eigenvec.txt"
    shell:
        "Rscript 'prepare/prepareKinship.R' --kinship={input.kinship}"
