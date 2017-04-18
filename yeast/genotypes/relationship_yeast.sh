###########################################################
###                                                     ###
###     Prune genotypes for SNPs in ld via plink 1.9    ###
###     * use different pruning parameters and for      ###
###       each set                                      ###
###         * compute pca                               ###
###         * compute genetic relationship              ###
###                                                     ###
###                                                     ###
###         by Hannah Meyer                             ###
###########################################################

dir=~/LiMMBo/yeast/inputdata
geno=BYxRM

### 1. Make binary file from ped/map
plink --file $dir/$geno --make-bed --out $dir/$geno
mv $dir/$geno.log $dir/log/$geno.makebed.log

# for different kb window sizes, do:
for kb in 3 5 8 10 12 15 20 30 50 100; do

    ### 2.  create ld-reduced , called genotype files
    plink --bfile $dir/$geno --indep-pairwise ${kb}kb 5 0.2 --out $dir/$geno.${kb}kb
    mv $dir/$geno.${kb}kb.log $dir/log/$geno.${kb}kb.prune.log

    plink --bfile $dir/$geno --extract $dir/$geno.${kb}kb.prune.in --make-bed --out $dir/$geno.${kb}kb.pruned
    mv $dir/$geno.${kb}kb.pruned.log $dir/log/$geno.${kb}kb.pruned.log

    ### 3. run pca
    plink --bfile $dir/$geno --extract $dir/$geno.prune.in --pca --out $dir/$geno.pca
    mv $dir/$geno.pca.log $dir/log/$geno.pca.log

    ### 4. compute genetic relationship matrix (grm)
    plink --bfile $dir/$geno --extract $dir/$geno.prune.in --make-rel square gz --out $dir/$geno.grm
    mv $dir/$geno.grm.log $dir/log/$geno.grm.log

    ### 5. format grm and depict as heatmap via population_structure.R
    Rscript ~/GWAS/analysis/GWAS/population_structure.R --vanilla  --default-packages=R.utils --data=$dir/$geno

    ### 6. compute LD and tag SNPs
    plink --bfile $dir/$geno --r2 with-freqs --show-tags all --tag-kb ${kb}kb --out $dir/$geno.${kb}kb
    mv $dir/$geno.${kb}kb.log $dir/log/$geno.computeLD${kb}kb.log
done
