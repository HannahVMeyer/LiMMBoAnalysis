#######################################################
###                                                 ###
###     Script to format .bed/.bim.fam to           ###
###                      .gen/.sample  to           ###
###                      limix_format               ###
###                                                 ###
###         by Hannah Meyer                         ###
#######################################################

dir=~/LiMMBo/yeast/inputdata
geno=BYxRM


# 1. Genotype
## a Create .gen/.sample files 
plink --bfile  $dir/$geno --recode oxford --out $dir/$geno
mv $dir/$geno.log $dir/log/$geno.genformat.log

## b. Reformat (to format needed as inpute for limix_format.py): chr22:16052684-rs139918843-A-C,A,C,2,1,1,2...
sed 's/_/ /g' $dir/$geno.gen | awk '{printf "chr"$1":"$4"-NA-"$5"-"$6","$8","$9; for(i=12;i<=NF;i+=3) printf ","$i; printf "\n"}' |gzip -f  > $dir/$geno.gen.gz

## c. Reformat for limix
~/GWAS/analysis/GWAS/limix_format_data.py -m genotype -f $dir/$geno.gen.gz -o $dir/$geno.geno.limix.format -s $dir/$geno.sample

# 2. Phenotype
~/GWAS/analysis/GWAS/limix_format_data.py -m phenotype -f $dir/${geno}_pheno_format.txt -o $dir/$geno.pheno.limix.format 

