#!/bin/bash 

##NJ method to build phylogenetic tree
##step1 data preparation--filter LD
plink --vcf  all_snp.recode.vcf  --indep-pairwise 50 5 0.5 --out tmp.ld   --allow-extra-chr --set-missing-var-ids @:#  
plink --vcf  all_snp.recode.vcf  --make-bed --extract tmp.ld.prune.in  --out all.LDfilter --recode vcf-iid  --keep-allele-order  --allow-extra-chr --set-missing-var-ids @:#
##step2 building phylogenetic tree
run_pipeline.pl -Xms1G -Xmx5G -importGuess LD_filter.vcf -ExportPlugin -saveAs sequences.phy -format Phylip_Inter
echo -e "sequences.phy\nY" > dnadist.cfg
dnadist < dnadist.cfg  >dnadist.log
mv outfile infile.dist
echo -e "infile.dist\nY"  > neighbor.cfg
neighbor  <  neighbor.cfg  >nj.log
less infile.dist | tr '\n' '|'| sed 's/| / /g' | tr '|' '\n' >infile.dist.table
less outtree | tr '\n' ' '|sed 's/ //g' > outtree.nwk

##PCA
plink --vcf  LD_filter.vcf   --pca 10 --out  PCA_out   --allow-extra-chr --set-missing-var-ids @:#

##Structure
plink --vcf  LD_filter.vcf  --make-bed --out all  --allow-extra-chr --keep-allele-order --set-missing-var-ids @:#
seq 2 10 | awk '{print "admixture --cv -j2 all.bed "$1" 1>admix."$1".log 2>&1"}' > admixture.sh
sh admixture.sh

##LD
/public/home/chcg/software/PopLDdecay/PopLDdecay -InVCF all_snp.recode.vcf -SubPop pop.txt -MaxDist 1 -OutStat pop.1k.stat ##pop.txt and HH.txt are sample-ids of different pop
ls pop.*.stat.gz |awk -F"." '{ print $0"\t"$2 }' > ld_stat.list
perl /public/home/chcg/software/PopLDdecay/bin/Plot_MultiPop.pl -inList ld_stat.list -output ld_stat.multi

##FST_Pi
python /public/home/chcg/software/genomics_general/VCF_processing/parseVCF.py -i all_snp.recode.vcf | gzip > subset.geno.gz
python /public/home/chcg/software/genomics_general/popgenWindows.py -g subset.geno.gz -o subset.Fst.Dxy.pi.csv.gz -f phased -w 10000  -s 10000 -p Cgi -p Can -p Csi -p Car --popsFile pop.info 
