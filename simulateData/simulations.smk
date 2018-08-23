# snakemake -s snake_simulations.smk --jobs 5000 --latency-wait 30 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources} -M {cluster.memory}  -o {cluster.output} -e {cluster.error}' --keep-going --rerun-incomplete
import numpy as np

configfile: "config/config_simulations.yaml"

rule all:
    input:
        expand("{genodir}/N10000/10000G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_kinship.csv",
            genodir=config["genodir"],
            A=config["paramGeno"]["nAnc"],
            B=config["paramGeno"]["nBlocks"],
            seed=config["paramGeno"]["seedGeno"],
            type=config["paramGeno"]["type"]),
        expand("{genodir}/N{N}/{N}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_chr{chr}_maf{maf}_gemma.csv",
            N=['1000', '10000'],
            genodir=config["genodir"],
            A=config["paramGeno"]["nAnc"],
            B=config["paramGeno"]["nBlocks"],
            seed=config["paramGeno"]["seedGeno"],
            type=config["paramGeno"]["type"],
            maf=config["paramGeno"]["maf"],
            chr=range(1,23)),
        expand("{phenodir}/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Ysim_plink.txt",
            phenodir=config["phenodir"],
            N=config["N"],
            seed=range(1,101),
            a=["{:.1f}".format(x) for x in np.append(0.05, np.arange(start=0.1,stop=1, step=0.1))],
            P=['50', '100', '1000'],
            h2=['0.3'],
            S=config["S"],
            model="modelnoiseFixedAndBggeneticFixedAndBg"),
        expand("{phenodir}/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Ysim_plink.txt",
            phenodir=config["phenodir"],
            N=1000,
            seed=range(101,201),
            a = "0.6",
            P=['10', '20', '30', '40', '50', '70', '100', '150', '200', '500', '1000'],
            h2=['0.3'],
            S=config["S"],
            model="modelnoiseFixedAndBggeneticFixedAndBg"),
        expand("{phenodir}/Calibration/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Ysim_plink.txt",
            phenodir=config["phenodir"],
            N=1000,
            seed=range(201,301),
            P=[10, 30, 100, 300],
            h2=['0.3'],
            S="0",
            model="modelnoiseFixedAndBggeneticBgOnly"),
        expand("{phenodir}/Scalability/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Ysim_plink.txt",
            phenodir=config["phenodir"],
            N=[1000, 2000, 5000, 10000],
            seed=range(201,301),
            P=100,
            h2=['0.3'],
            S="0",
            model="modelnoiseFixedAndBggeneticBgOnly"),

rule simulateGenotypes:
    input:
        Genomes=lambda wildcards: config["Genomes"]
    params:
        calc_covar="True",
        maf=config["paramGeno"]["maf"]
    output:
        "{dir}/genotypes/relatedEU_nopopstructure/N{N}/{N}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_kinship.csv"
    log:
        "{dir}/genotypes/relatedEU_nopopstructure/N{N}/log/{N}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_kinship.log"
    shell:
        "(python genotypes/simulateGenotypes.py \
            --nSamples {wildcards.N} --nAncestors {wildcards.A} --chroms 1-22 \
            --seed 256 \
            --out_dir {wildcards.dir}/relatedEU_nopopstructure/N{wildcards.N} \
            --type='noPop' --maf {params.maf} --csv) 2> {log}"

rule plotSimulateGenotypes:
    input:
        "{dir}/genotypes/relatedEU_nopopstructure/N{N}/{N}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_kinship.csv"
    output:
        "{dir}/genotypes/relatedEU_nopopstructure/N{N}/simulatedCovarianceMatrices_kinship.png"
    log:
        "{dir}/genotypes/relatedEU_nopopstructure/N{N}/simulatedCovarianceMatrices_kinship.log"
    shell:
        "(genotypes/plotSimulateGenotypes.R) 2> {log}"

rule formatGenotyoes:
    input:
        csv="{dir}/genotypes/relatedEU_nopopstructure/N{N}/{N}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_chr{chr}_maf{maf}.csv"
    output:
        "{dir}/genotypes/relatedEU_nopopstructure/N{N}/{N}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_chr{chr}_maf{maf}_gemma.csv"
    shell:
        "awk -F, '{{$1=$1\",A,T\";}}1' OFS=, {input.csv} > {output}"

rule simulateTraitsAffected:
    input:
        kinship=expand("{{dir}}/genotypes/relatedEU_nopopstructure/N{N}/{N}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_kinship_norm.csv",
            N=config["N"], A=config["paramGeno"]["nAnc"],
            B=config["paramGeno"]["nBlocks"],
            type=config["paramGeno"]["type"],
            seed=config["paramGeno"]["seedGeno"])
    params:
        genoPrefix=expand("genotypes/relatedEU_nopopstructure/N{N}/{N}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_",
            N=config["N"], A=config["paramGeno"]["nAnc"],
            B=config["paramGeno"]["nBlocks"],
            seed=config["paramGeno"]["seedGeno"],
            type=config["paramGeno"]["type"]),
        genoSuffix="_maf0.02.csv",
        h2s="0.05"
    output:
        "{dir}/phenotypes/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Ysim_plink.txt"
    log:
        "{dir}/phenotypes/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/log/sim.log"
    shell:
        "Rscript -e 'PhenotypeSimulator::simulatePhenotypes()' \
            --args \
            --NrSamples={wildcards.N} --NrPhenotypes={wildcards.P} \
            --cNrSNP={wildcards.S} \
            --genoFilePrefix={wildcards.dir}/{params.genoPrefix} \
            --genoFileSuffix={params.genoSuffix} \
            --NrChrCausal=5 \
            --kinshipFile={input.kinship} \
            --pTraitsAffectedGenetics={wildcards.a} \
            --seed={wildcards.seed} \
            --genVar={wildcards.h2} --h2s={params.h2s} \
            --delta=0.4 --phi=0.6 \
            --showProgress \
            --subdirectory="" \
            --directory={wildcards.dir}/phenotypes/TraitsAffected{wildcards.a}/Traits{wildcards.P}_samples{wildcards.N}_NrSNP{wildcards.S}_Cg{wildcards.h2}_modelnoiseFixedAndBggeneticFixedAndBg/seed{wildcards.seed} \
            --saveTable \
            --savePLINK \
            --saveSNPTEST \
            --saveGEMMA"

rule simulateCalibration:
    input:
        kinship=expand("{{dir}}/genotypes/relatedEU_nopopstructure/N{N}/{N}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_kinship_norm.csv",
            N=config["N"], A=config["paramGeno"]["nAnc"],
            B=config["paramGeno"]["nBlocks"],
            type=config["paramGeno"]["type"],
            seed=config["paramGeno"]["seedGeno"])
    params:
        h2s="0"
    output:
        "{dir}/phenotypes/Calibration/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Ysim_plink.txt"
    log:
        "{dir}/phenotypes/Calibration/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/log/sim.log"
    shell:
        "Rscript -e 'PhenotypeSimulator::simulatePhenotypes()' \
            --args \
            --NrSamples={wildcards.N} --NrPhenotypes={wildcards.P} \
            --cNrSNP={wildcards.S} \
            --kinshipFile={input.kinship} \
            --seed={wildcards.seed} \
            --genVar={wildcards.h2} --h2s={params.h2s} \
            --delta=0.4 --phi=0.6 \
            --showProgress \
            --subdirectory="" \
            --directory={wildcards.dir}/phenotypes/Calibration/Traits{wildcards.P}_samples{wildcards.N}_NrSNP{wildcards.S}_Cg{wildcards.h2}_{wildcards.model}/seed{wildcards.seed} \
            --saveTable \
            --savePLINK \
            --saveSNPTEST \
            --saveGEMMA"

rule simulateScalabilitySamples:
    input:
        kinship=expand("{{dir}}/genotypes/relatedEU_nopopstructure/N{N}/{N}G_nAnc{A}_nBlocks{B}_seed{seed}_{type}_kinship_norm.csv",
            N=10000, A=config["paramGeno"]["nAnc"],
            B=config["paramGeno"]["nBlocks"],
            type=config["paramGeno"]["type"],
            seed=config["paramGeno"]["seedGeno"])
    params:
        h2s="0"
    output:
        "{dir}/phenotypes/Scalability/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Ysim_plink.txt",
    log:
        "{dir}/phenotypes/Scalability/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/log/sim.log",
    shell:
        "Rscript -e 'PhenotypeSimulator::simulatePhenotypes()' \
            --args \
            --NrSamples={wildcards.N} --NrPhenotypes={wildcards.P} \
            --cNrSNP={wildcards.S} \
            --kinshipFile={input.kinship} \
            --seed={wildcards.seed} \
            --genVar={wildcards.h2} --h2s={params.h2s} \
            --delta=0.4 --phi=0.6 \
            --showProgress \
            --subdirectory="" \
            --directory={wildcards.dir}/phenotypes/Scalability/Traits{wildcards.P}_samples{wildcards.N}_NrSNP{wildcards.S}_Cg{wildcards.h2}_{wildcards.model}/seed{wildcards.seed} \
            --saveTable \
            --savePLINK \
            --saveSNPTEST \
            --saveGEMMA"
