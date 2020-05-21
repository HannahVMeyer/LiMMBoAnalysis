# snakemake -s power.smk --jobs 5000 --latency-wait 30 --cluster-config cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources} -M {cluster.memory}  -o {cluster.output} -e {cluster.error}' --keep-going --rerun-incomplete
import numpy as np

configfile: "config_power.yaml"

rule all:
    input:
        expand("{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/{type}",
            dir=config["dir"],
            N='1000',
            seed=range(1,101),
            a=np.append(["0.05"], ["{:.1f}".format(x) for x in np.append(np.arange(start=0.1,stop=0.6, step=0.1), np.arange(start=0.7,stop=1, step=0.1))]),
            P=['50', '100', '1000'],
            h2=['0.3'],
            S=config["S"],
            model="modelnoiseFixedAndBggeneticFixedAndBg",
            type=['sbat.C.mpmm.txt', 'Cg_gemma.csv']),
        expand("{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/{type}",
            dir=config["dir"],
            N=['5000', '10000']
            seed=range(101,201),
            a = "0.6",
            P=['10', '20', '30', '40', '50', '70', '100', '150', '200', '500', '1000'],
            h2=['0.3'],
            S=config["S"],
            model="modelnoiseFixedAndBggeneticFixedAndBg",
            type=['sbat.C.mpmm.txt', 'Cg_gemma.csv'])



rule power_limmbo:
    input:
        eigenvec=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/kinship_eigenvec_limmbo.csv",
            dir=config["simulatedir"]),
        eigenval=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/kinship_eigenvalue_limmbo.csv",
            dir=config["simulatedir"]),
        pheno=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Ysim_reg_limmbo.csv",
            dir=config["simulatedir"])
    params:
        sp=lambda wildcards: '5' if wildcards.P == 10 else '10',
        min=config["params"]["min"],
        tr=config["params"]["tr"],
        cpus=config["params"]["cpu"]
    output:
        "{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Cg_fit.csv",
        "{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Cn_fit.csv"
    log:
        "{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/log/limmbo.log"
    conda:
        "envs/limmbo1.0.0-dev.yaml"
    shell:
        "(runVarianceEstimation \
            -o {wildcards.dir}/vd/TraitsAffected{a}/Traits{wildcards.P}_samples{wildcards.N}_NrSNP{wildcards.S}_Cg{wildcards.h2}_{wildcards.model}/seed{wildcards.seed} \
            -sp {params.sp} -seed {wildcards.seed} \
            -p {input.pheno} \
            -tr {params.tr} \
            -eval_k {input.eigenval} \
            -evec_k {input.eigenvec} \
            --minCooccurrence {params.min} \
            -cpus {params.cpus} --limmbo -v -t) 2> {log}"

rule power_limix:
    input:
        eigenvec=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/kinship_eigenvec_limmbo.csv",
            dir=config["simulatedir"]),
        eigenval=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/kinship_eigenvalue_limmbo.csv",
            dir=config["simulatedir"]),
        pheno=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Ysim_reg_limmbo.csv",
            dir=config["simulatedir"])
    params:
        tr=config["params"]["tr"],
    output:
        "{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Cg_REML.csv",
        "{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Cn_REML.csv",
    log:
        "{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/log/limix.log"
    conda:
        "envs/limmbo1.0.0-dev.yaml"
    shell:
        "(runVarianceEstimation \
            -o {wildcards.dir}/vd/TraitsAffected{a}/Traits{wildcards.P}_samples{wildcards.N}_NrSNP{wildcards.S}_Cg{wildcards.h2}_{wildcards.model}/seed{wildcards.seed} \
            -seed {wildcards.seed} \
            -p {input.pheno} \
            -eval_k {input.eigenval} \
            -evec_k {input.eigenvec} \
            -tr {params.tr} \
            --reml -v -t) 2> {log}"

rule power_sbat:
    input:
        eigenvec=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/kinship_eigenvec.csv",
            dir=config["simulatedir"]),
        eigenval=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/kinship_eigenval.csv",
            dir=config["simulatedir"]),
        pheno=expand("{dir}/phenotypes/Calibration/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Ysim_reg_sbat.txt",
            dir=config["simulatedir"])
    output:
        C="{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/sbat.C.mpmm.txt",
        D="{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/sbat.D.mpmm.txt",
        sbatlog="{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/sbat.log",
        time="{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/VD_time_sbat.csv"
    log:
        "{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/log/sbat.log"
    shell:
        "(sbat_fs --covariances-only \
              -p {input.pheno} \
              -ev {input.eigenvec} {input.eigenval} \
              -o {wildcards.dir}/vd/TraitsAffected{a}/Traits{wildcards.P}_samples{wildcards.N}_NrSNP{wildcards.S}_Cg{wildcards.h2}_{wildcards.model}/seed{wildcards.seed}/sbat) 2> {log}"

rule power_gemma:
    input:
        eigenvec=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/kinship_eigenvec.csv",
            dir=config["simulatedir"]),
        eigenval=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/kinship_eigenval.csv",
            dir=config["simulatedir"]),
        pheno=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Ysim_reg_gemma.txt",
            dir=config["simulatedir"])
    output:
        gemmalog="{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/gemma.log",
        time="{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/VD_time_gemma.csv",
        ve="{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Cn_gemma.csv",
        vg="{dir}/vd/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/Cg_gemma.csv"
    params:
        num=lambda wildcards: ' '.join(map(str, range(1,(int(wildcards.P)+1))))
    shell:
        """
        (gemma -p {input.pheno} \
           -g {input.geno} \
           -d {input.eigenval} \
           -u {input.eigenvec} \
           -lmm 2 \[
           -n {params.num} \
           -o {wildcards.dir}/vd/TraitsAffected{wildcards.a}/Traits{wildcards.P}_samples{wildcards.N}_NrSNP{wildcards.S}_Cg{wildcards.h2}_{wildcards.model}/seed{wildcards.seed}) 2> {output.gemmalog}
        GemmaVD-parser.py -f {output.log} -o {wildcards.dir}/vd/TraitsAffected{wildcards.a}/Traits{wildcards.P}_samples{wildcards.N}_NrSNP{wildcards.S}_Cg{wildcards.h2}_{wildcards.model}/seed{wildcards.seed}
        """

rule association_limmbo:
    input:
        eigenvec=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/kinship_eigenvec_limmbo.csv",
            dir=config["simulatedir"]),
        eigenval=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/kinship_eigenvalue_limmbo.csv",
            dir=config["simulatedir"]),
        pheno=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Ysim_reg_limmbo.csv",
            dir=config["simulatedir"]),
        geno=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/SNP_NrSNP20_limmbo.csv",
            dir=config["simulatedir"]),
        cg=expand("{dir}/vd/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Cg_fit.csv",
            dir=config["dir"]),
        cn=expand("{dir}/vd/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Cn_fit.csv",
            dir=config["dir"])
    params:
        tr=config["params"]["tr"]
    output:
        "{dir}/association/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/lmm_mt_pvalues_limmbo.csv"
    log:
        "{dir}/association/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/log/lmm_mt_limmbo.log"
    conda:
        "envs/limmbo1.0.0-dev.yaml"
    shell:
        "(runAssociation \
            -o {wildcards.dir}/association/TraitsAffected{a}/Traits{wildcards.P}_samples{wildcards.N}_NrSNP{wildcards.S}_Cg{wildcards.h2}_{wildcards.model}/seed{wildcards.seed} \
            -p {input.pheno} \
            -cg {input.cg} \
            -cn {input.cn} \
            -eval_k {input.eigenval} \
            -evec_k {input.eigenvec} \
            -tr {params.tr} \
            -v -t) 2> {log}"

rule association_limix:
    input:
        eigenvec=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/kinship_eigenvec_limmbo.csv",
            dir=config["simulatedir"]),
        eigenval=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/kinship_eigenvalue_limmbo.csv",
            dir=config["simulatedir"]),
        pheno=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Ysim_reg_limmbo.csv",
            dir=config["simulatedir"]),
        geno=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/SNP_NrSNP20_limmbo.csv",
            dir=config["simulatedir"]),
        cg=expand("{dir}/vd/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Cg_REML.csv",
            dir=config["dir"]),
        cn=expand("{dir}/vd/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Cn_REML.csv",
            dir=config["dir"])
    params:
        tr=config["params"]["tr"],
    output:
        "{dir}/association/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/lmm_mt_pvalues_limix.csv"
    log:
        "{dir}/association/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/log/lmm_mt_limix.log"
    conda:
        "envs/limmbo1.0.0-dev.yaml"
    shell:
        "(runAssociation \
            -o {wildcards.dir}/vd/TraitsAffected{a}/Traits{wildcards.P}_samples{wildcards.N}_NrSNP{wildcards.S}_Cg{wildcards.h2}_{wildcards.model}/seed{wildcards.seed} \
            -seed {wildcards.seed} \
            -p {input.pheno} \
            -g {input.geno} \
            -cg {input.cg} \
            -cn {input.cn} \
            -eval_k {input.eigenval} \
            -evec_k {input.eigenvec} \
            -tr {params.tr} \
            -v -t) 2> {log}"

rule association_sbat:
    input:
        eigenvec=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/kinship_eigenvec.csv",
            dir=config["simulatedir"]),
        eigenval=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/kinship_eigenval.csv",
            dir=config["simulatedir"]),
        pheno=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Ysim_reg_sbat.txt",
            dir=config["simulatedir"])
        geno=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/SNP_NrSNP20_snptest.gen",
            dir=config["simulatedir"]),
        cg=expand("{dir}/vd/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/sbat.C.mpmm.txt",
            dir=config["simulatedir"]),
        cn=expand("{dir}/vd/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/sbat.C.mpmm.txt",
            dir=config["simulatedir"])
    output:
        "{dir}/association/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/sbat.MLMM.summary.txt"
    log:
        "{dir}/association/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/sbat.MLMM.summary.log"
    shell:
        "(sbat_fs --full-model \
              -p {input.pheno} \
              -g {input.geno} \
              -C {input.cg} \
              -D {input.cn} \
              -o {wildcards.dir}/vd/TraitsAffected{a}/Traits{wildcards.P}_samples{wildcards.N}_NrSNP{wildcards.S}_Cg{wildcards.h2}_{wildcards.model}/seed{wildcards.seed}/sbat) 2> {log}"

rule association_gemma:
    input:
        pheno=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/Ysim_reg_gemma.txt",
            dir=config["simulatedir"]),
        geno=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/SNP_NrSNP20_gemma.txt",
            dir=config["simulatedir"]),
        eigenvec=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/kinship_eigenvec.csv",
            dir=config["simulatedir"]),
        eigenval=expand("{dir}/phenotypes/TraitsAffected{{a}}/Traits{{P}}_samples{{N}}_NrSNP{{S}}_Cg{{h2}}_{{model}}/seed{{seed}}/kinship_eigenval.csv",
            dir=config["simulatedir"]),
    output:
        "{dir}/association/TraitsAffected{a}/Traits{P}_samples{N}_NrSNP{S}_Cg{h2}_{model}/seed{seed}/gemma.assoc.txt"
    params:
        num=lambda wildcards: ' '.join(map(str, range(1,(int(wildcards.P)+1))))
    shell:
        """
        (gemma -p {input.pheno} \
           -g {input.geno} \
           -d {input.eigenval} \
           -u {input.eigenvec} \
           -lmm 2 \
           -n {params.num} \
           -o {wildcards.dir}/vd/TraitsAffected{wildcards.a}/Traits{wildcards.P}_samples{wildcards.N}_NrSNP{wildcards.S}_Cg{wildcards.h2}_{wildcards.model}/seed{wildcards.seed}/gemma) 2> {log}
        """
