#!/bin/bash 
##This script is used to call SNPs based on 95 oyster resequencing samples
###step1 Filtering and quality control of raw data
fastp -i  ${i}_1_clean.fq.gz -I ${i}_2_clean.fq.gz -o ${i}_1_out.fq.gz -O ${i}_2_out.fq.gz -h out_R.html -j out_R.json
###step2 mapping and sort 
bwa mem -t 20 -R '@RG\tID:'$i'\tSM:'$i'\tPL:illumina'  ref_genome.fasta  ${i}_1_clean_val_1.fq.gz  ${i}_2_clean_val_2.fq.gz | samtools sort -@ 20 -m 1G  -o  ${i}.sort.bam -
###step3 remove PCR repetation
picard -Xmx4g MarkDuplicates I=bei_${i}.sort.bam O=bei_${i}.sort.rmdup.bam CREATE_INDEX=true REMOVE_DUPLICATES=true   M= bei_${i}.sort.markdup_metrics.txt
###step4 haplotyecaller
gatk --java-options "-Xmx5g -Djava.io.tmpdir=./tmp" HaplotypeCaller -R  refer_genome.fa -I bei_${i}.sort.rmdup.bam  -ERC GVCF -O  bei_${i}.g.vcf 1>bei_${i}.HC.log   2>&1  #1.need create a tmp file to save temporary documents
###step5 merge all gvcf files
ls  bei_*.g.vcf > ./gvcf.list
gatk  --java-options "-Xmx4g -Djava.io.tmpdir=./tmp"   CombineGVCFs -R  refer_gemoe.fa  -V ./gvcf.list  -O ./all.merge.g.vcf
###step6 population variation calling 
gatk  --java-options "-Xmx100g -Djava.io.tmpdir=./tmp"   GenotypeGVCFs -R refer_genome.fa --variant  all.merge.g.vcf -O ./all.merge_raw.vcf
###step7 extract SNP
gatk  --java-options "-Xmx100g  -Djava.io.tmpdir=./tmp"  SelectVariants  -R refer_genome.fa -V all.merge_raw.vcf  --select-type SNP -O ./all.raw.snp.vcf
###step8 hard filter
gatk  --java-options "-Xmx100g -Djava.io.tmpdir=./tmp"  VariantFiltration -R refer_genome.fa  -V  ./all.raw.snp.vcf  --filter-expression "QUAL <50 || QD < 2.0 ||  DP < 3 || DP > 100  || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name 'SNP_filter' -O ./all.filter.snp.vcf
gatk --java-options "-Xmx100g -Djava.io.tmpdir=./tmp" SelectVariants -R .refer_genome.fa -V all.filter.snp.vcf --exclude-filtered -O all.filtered.snp.vcf
###step9 soft filter
vcftools --vcf  all.filter.snp.vcf  --recode-INFO-all   --max-alleles 2   --min-alleles 2  --max-missing 0.8 --out all_snp --recode --remove-filtered-all  # other parameters --maf 0.05 --minDP 3 --min-meanDP 3 --minQ 200 
