#### ## redundant_gene
cp /public1/node3_liushk/tmp/03.braker_result/Cskiamea*pep.fa .
conda activate genome

cd-hit -i CskiameaA.pep.fa  -o HAPA_90 -c 0.9 -n 5 -g 1 -G 0 -aS 0.8  -d 0 -p 1 -T 16 -M 0
cd-hit -i CskiameaB.pep.fa  -o HAPB_90 -c 0.9 -n 5 -g 1 -G 0 -aS 0.8  -d 0 -p 1 -T 16 -M 0

clstr2txt.pl HAPA_90.clstr  > hapA.clstr.txt
clstr2txt.pl HAPB_90.clstr  > hapB.clstr.txt

cut -f 2 hapA.clstr.txt | sort | uniq -c | awk '{if($1>1)print $2}' > cluts.tmp
while read i ; do awk '{if($2=="'$i'") print $0}' hapA.clstr.txt >> hapA.tmp ; done < cluts.tmp
cut -f 1 hapA.tmp  > 01.hapA.redundancy.gene

cut -f 2 hapB.clstr.txt | sort | uniq -c | awk '{if($1>1)print $2}' > cluts.tmp
while read i ; do awk '{if($2=="'$i'") print $0}' hapB.clstr.txt >> hapB.tmp ; done < cluts.tmp
cut -f 1 hapB.tmp  > 01.hapB.redundancy.gene

bedtools coverage -a 01.hapA.redundancy.gene.gff -b 00.data/03.noalign.bed  | awk '{if($13==1) print $0}' | wc -l
bedtools coverage  -a 00.all.gene.gff -b 00.data/03.syn.bed  | awk '{if($13==1) print $0}' | wc -l

#### same pipeline for other species' genome


## plot in R
# 创建数据框
data <- data.frame(
    group = c("HDS-csi", "SYNAL-csi", "SVAL-csi",
              "HDS-ase", "SYNAL-ase", "SVAL-ase",
              "HDS-pfu", "SYNAL-pfu", "SVAL-pfu",
              "HDS-mva", "SYNAL-mva", "SVAL-mva"),
    genome_dup = c(4326, 4326, 4326, 8363, 8363, 8363,
                   3967, 3967, 3967, 5001, 5001, 5001),
    genome_total = c(30052, 30052, 30052, 33651, 33651, 33651,
                     32938, 32938, 32938, 29936, 29936, 29936),
    region_total = c(11131, 4266, 570, 16657, 1691, 774,
                     11686, 1845, 253, 9825, 3843, 492),
    region_dup = c(3035, 379, 247, 5967, 210, 425,
                   2669, 165, 106, 3543, 461, 186),
    stringsAsFactors = FALSE
)

# 检查数据类型
str(data)

# 计算 p 值并添加到新列
data$p_value <- apply(data, 1, function(x) {
    # 提取每行的数值
    genome_dup <- as.numeric(x["genome_dup"])
    genome_total <- as.numeric(x["genome_total"])
    region_total <- as.numeric(x["region_total"])
    region_dup <- as.numeric(x["region_dup"])

    # 计算特定区域和非特定区域的重复和非重复基因数
    overlap_genes <- region_dup
    non_overlap_genes <- region_total - region_dup
    genome_dup_genes <- genome_dup - overlap_genes
    genome_non_dup_genes <- genome_total - genome_dup - non_overlap_genes

    # 构建 2x2 列联表
    contingency_table <- matrix(c(overlap_genes, non_overlap_genes,
                                  genome_dup_genes, genome_non_dup_genes),
                                nrow = 2,
                                byrow = TRUE)

    # 检查表格是否合理
    if(any(contingency_table < 0)) {
        warning(paste("负数在列联表中发现，检查数据行:", x["group"]))
        return(NA)
    }

    # 执行 Fisher 精确检验，选择 'greater' 假设（检测富集）
    test <- fisher.test(contingency_table, alternative = "greater")

    return(test$p.value)
})

# 查看结果
print(data)
data$padj <- p.adjust(data$p_value, method = "BH")

data$Rich.factor <- data$region_dup/data$genome_dup

library(ggplot2)
library(forcats)

A<- data
# 替代 0 值为 1e-10
A$padj <- ifelse(data$padj == 0, 1e-100, A$padj)

# 对 p 值和 padj 进行负对数转换
A$log_padj <- -log10(A$padj)

A$Description <- as.factor(A$group)
A$Description <- fct_inorder(A$Description)
sub(".*-(.*)", "\\1", A$group) -> A$group1
sub("-(.*)", "", A$group) ->  A$group2

A$group1 <- factor(A$group1,levels = c("csi", "pfu", "ase","mva"))
ggplot(A, aes(group1, group2)) +
    geom_point(aes(color=log_padj, size=Rich.factor))+theme_bw()+
    theme(panel.grid = element_blank(),
          axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
    scale_color_gradient(low='#6699CC',high='#CC3333')+
    labs(x=NULL,y=NULL)+guides(size=guide_legend(order=1)) -> p

p1 <- p + guides(colour = guide_colorbar(order = 1),size = guide_legend(order = 2))
p2<- p1 + theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=14,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))

pdf("plot2kegg.pdf", width=4, height=3)
p2 +coord_flip()
dev.off()

#############################################################################
### 00.trimmed and filtered
while read i ; do fastp -i $i'_1.fq.gz' -I $i'_2.fq.gz' -o $i'cleaned_1.fq.gz' -O $i'cleaned_2.fq.gz' -l 36 -q 20 -w 12 ; done < list

### 01.reads_mapping
STAR v2.7.9a was used with default parameters for the read mapping on the C.gigas reference genome with GTF file, by the Multi-sample 2-pass method.
1. Run 1st mapping pass for all samples with ”usual” parameters. Using annotations is recommended either a the genome generation step, or mapping step.
2. Run 2nd mapping pass for all samples , listing SJ.out.tab files from all samples in
--sjdbFileChrStartEnd /path/to/sj1.tab /path/to/sj2.tab ....

# 01./1/ star 1-pass index
cat HAP1.fasta HAP2.fasta > genome.fa
cat /public1/node3_liushk/tmp/03.braker_result/CskiameaA.longest.gff  /public1/node3_liushk/tmp/03.braker_result/CskiameaB.longest.gff > Cskiamea.d.gff
agat_convert_sp_gff2gtf.pl --gff Cskiamea.d.gff -o Cskiamea.d.gtf #genomic
STAR --runMode genomeGenerate --runThreadN 16 --genomeDir 00.star_index/ --genomeFastaFiles genome.fa --sjdbGTFfile Cskiamea.d.gtf --sjdbOverhang 149 #genome

# 01./2/ 1-pass mapping with indexed genome
while read i ; do STAR --genomeDir 00.star_index/ --runThreadN 20  --readFilesIn '01.fastp/'$i'_1.cleaned.fq.gz'  '01.fastp/'$i'_2.cleaned.fq.gz'  --readFilesCommand zcat --outSAMunmapped Within --outFileNamePrefix '02.1-pass/'$i  > '02.1-pass/log2map.'$i 2>&1  ; done < 01.fastp/list
cat *.tab | awk '($5 > 0 && $7 > 2 && $6==0)' | cut -f1-6 | sort | uniq > SJ.filtered.tab

# 2-pass mapping with SJ.out.tab files
while read i ; do STAR --genomeDir 00.star_index/ --runThreadN 30 --sjdbFileChrStartEnd SJ.filtered.tab --readFilesCommand zcat \
--readFilesIn '01.fastp/'$i'_1.cleaned.fq.gz' '01.fastp/'$i'_2.cleaned.fq.gz' \
--outSAMmapqUnique 255 --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 30 \
--outFileNamePrefix '03.2-pass/'$i  > '03.2-pass/log2map.'$i 2>&1 ; done <  01.fastp/list


**注意，该部分数据用于 个体的等位基因表达和ASE分析**

### /public1/node3_shicy/PRJ02_Csikamea/00.ASE
cat CskiameaA.longest.gff CskiameaB.longest.gff > Cskiamea.d.gff
stringtie -p 20 -G Cskiamea.d.gff  -e -B -o 04.string/TD4/transcripts.gtf -A 04.string/TD4/gene_abundances.tsv 03.2-pass/TD4Aligned.sortedByCoord.out.bam
stringtie -p 20 -G Cskiamea.d.gff  -e -B -o 04.string/TG4/transcripts.gtf -A 04.string/TG4/gene_abundances.tsv 03.2-pass/TG4Aligned.sortedByCoord.out.bam
stringtie -p 20 -G Cskiamea.d.gff  -e -B -o 04.string/TL4/transcripts.gtf -A 04.string/TL4/gene_abundances.tsv 03.2-pass/TL4Aligned.sortedByCoord.out.bam
stringtie_expression_matrix.pl --expression_metric=TPM --result_dirs='TD4,TG4,TL4' --transcript_matrix_file=transcript_tpm_all_samples.tsv --gene_matrix_file=gene_tpm_all_samples.tsv

## /public1/node3_shicy/PRJ02_Csikamea/00.ASE/05.expr2region
bedtools coverage  -a /public1/node3_shicy/PRJ02_Csikamea/12.redundant_gene/00.all.gene.gff -b /public1/node3_shicy/PRJ02_Csikamea/12.redundant_gene/00.data/03.noalign.bed | awk '{if($13==1) print $0}' > 01.notalign.gene.bed
bedtools coverage  -a /public1/node3_shicy/PRJ02_Csikamea/12.redundant_gene/00.all.gene.gff -b /public1/node3_shicy/PRJ02_Csikamea/12.redundant_gene/00.data/03.sv.bed | awk '{if($13==1) print $0}' > 01.sv.gene.bed
bedtools coverage  -a /public1/node3_shicy/PRJ02_Csikamea/12.redundant_gene/00.all.gene.gff -b /public1/node3_shicy/PRJ02_Csikamea/12.redundant_gene/00.data/03.syn.bed | awk '{if($13==1) print $0}' > 01.sy.gene.bed

awk '{if($1~"A")print $1"\t"$9"\tHAPA";else print $1"\t"$9"\tHAPB"}' 01.notalign.gene.bed | sed 's/Achr/chr/g' | sed 's/Bchr/chr/g' | sed 's/ID=//g' | awk '{print $0"\t""NAL"}' > 02.not.gene
awk '{if($1~"A")print $1"\t"$9"\tHAPA";else print $1"\t"$9"\tHAPB"}' 01.sv.gene.bed       | sed 's/Achr/chr/g' | sed 's/Bchr/chr/g' | sed 's/ID=//g' | awk '{print $0"\t""SV"}' >  02.SV.gene
awk '{if($1~"A")print $1"\t"$9"\tHAPA";else print $1"\t"$9"\tHAPB"}' 01.sy.gene.bed       | sed 's/Achr/chr/g' | sed 's/Bchr/chr/g' | sed 's/ID=//g' | awk '{print $0"\t""SY"}' >  02.SY.gene
awk '{if($1~"A")print $1"\t"$9"\tHAPA";else print $1"\t"$9"\tHAPB"}' /public1/node3_shicy/PRJ02_Csikamea/12.redundant_gene/00.all.gene.gff | sed 's/Achr/chr/g' | sed 's/Bchr/chr/g' | sed 's/ID=//g'  | awk '{print $0"\t""ALL"}' >  02.all.gene

cp /public1/node3_shicy/PRJ02_Csikamea/00.ASE/04.string/gene_tpm_all_samples.tsv 03.gene.tpm.tsv


##### in R  Fig2b

library(plyr)
library(ggthemes)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(scales)
library(gghalves)

###自定义颜色
mypal=pal_simpsons(alpha = .6)(9)
mypal[c(1,7)]
show_col(mypal)
show_col(mypal[c(1,7)])

###自定义主题
mytheme <- theme(axis.text.x=element_text(size=12),
                 axis.text.y=element_text(size=12),
                 axis.title=element_text(size = 13),
                 legend.text=element_text(size=12),
                 legend.title=element_text(size=12),
                 axis.line = element_line(linewidth=0.7),
                 panel.border = element_blank(),
                 panel.grid = element_blank())

read.table('gene_tpm_all_samples.tsv',header = 1) -> aa_all
read.table('00.hapa.HDA.gene') -> hds_id
colnames(hds_id) <- "Gene_ID"
merge(hds_id,aa_all,by="Gene_ID") -> hds_expr   

read.table('00.geneinSV.csi') -> hds_id
colnames(hds_id) <- "Gene_ID"
merge(hds_id,aa_all,by="Gene_ID") -> sv_expr   

read.table('00.geneinSYNAL.csi') -> hds_id
colnames(hds_id) <- "Gene_ID"
merge(hds_id,aa_all,by="Gene_ID") -> sy_expr   


TD.hds.df <- data.frame("Group"="hds.TD","value"=hds_expr$TD4)
TG.hds.df <- data.frame("Group"="hds.TG","value"=hds_expr$TG4)
TL.hds.df <- data.frame("Group"="hds.TL","value"=hds_expr$TL4)
TD.sy.df <- data.frame("Group"="sy.TD","value"=sy_expr$TD4)
TG.sy.df <- data.frame("Group"="sy.TG","value"=sy_expr$TG4)
TL.sy.df <- data.frame("Group"="sy.TL","value"=sy_expr$TL4)
TD.sv.df <- data.frame("Group"="sv.TD","value"=sv_expr$TD4)
TG.sv.df <- data.frame("Group"="sv.TG","value"=sv_expr$TG4)
TL.sv.df <- data.frame("Group"="sv.TL","value"=sv_expr$TL4)


rbind(TD.hds.df,TG.hds.df,TL.hds.df,TD.sy.df,TG.sy.df,TL.sy.df,TD.sv.df,TG.sv.df,TL.sv.df) -> df

source("geom_flat_violin.R")

p <- ggplot(df,aes(x=Group,y=log2(value+1), fill = Group, colour = Group))+
geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2)+
geom_boxplot(aes(x = Group, y = log2(value+1), fill = Group),outlier.shape = NA, alpha = .5, width = .1, colour = "black")
comparisons <- list(c("hds.TD", "sv.TD"),c("hds.TD", "sy.TD"),c("sv.TD", "sy.TD"),c("hds.TG", "sv.TG"),c("hds.TG", "sy.TG"),c("sv.TG", "sy.TG"),c("hds.TL", "sv.TL"),c("hds.TL", "sy.TL"),c("sv.TL", "sy.TL"))

pdf(file = "fig4b.pdf", width =6, height = 4)
p+stat_compare_means(comparisons=comparisons,
                            label = "p.signif",method = "wilcox.test",hide.ns = F) + theme(legend.position = "none") + mytheme
dev.off( )


#######################################################################################################################################


orthinder for new gene identify

awk -F "\t" '{if($9>0 && $10>0)print $0}'  Orthogroups.GeneCount.tsv | awk -F "\t" '{if($4>0 && $5>0 && $7>0 && $8>0 && $11>0)print $0}' | awk -F "\t" '{if($17>0 && $18>0)print $0}' |  awk -F "\t" '{if($20>0 && $13>0 && $15>0 && $14>0 && $3>0 && $6>0 && $21>0 ) print $0}' | awk -F "\t" '{if($22>0 && $23>0 && $12>0 ) print $0}' | awk -F "\t" '{if($2==0 && $19==0 && $16==0 ) print $0}' | cut -f 1 > 01.Bivalvia

awk -F "\t" '{if($9>0 && $10>0)print $0}'  Orthogroups.GeneCount.tsv | awk -F "\t" '{if($4>0 && $5>0 && $7>0 && $8>0 && $11>0)print $0}' | awk -F "\t" '{if($17>0 && $18>0)print $0}' |  awk -F "\t" '{if($20>0 && $13>0 && $15>0 && $14>0 && $3>0 && $6>0 && $21>0 ) print $0}' | awk -F "\t" '{if($22==0 && $23==0 && $12==0 ) print $0}' | awk -F "\t" '{if($2==0 && $19==0 && $16==0 ) print $0}' | cut -f 1 >  02.Pteriomorphia

awk -F "\t" '{if($9>0 && $10>0)print $0}'  Orthogroups.GeneCount.tsv | awk -F "\t" '{if($4>0 && $5>0 && $7>0 && $8>0 && $11>0)print $0}' | awk -F "\t" '{if($17>0 && $18>0)print $0}' |  awk -F "\t" '{if($20==0 && $13==0 && $15==0 && $14==0 && $3==0 && $6==0 && $21==0 ) print $0}' | awk -F "\t" '{if($22==0 && $23==0 && $12==0 ) print $0}' | awk -F "\t" '{if($2==0 && $19==0 && $16==0 ) print $0}' | cut -f 1 >  03.Osteridae

awk -F "\t" '{if($9>0 && $10>0)print $0}'  Orthogroups.GeneCount.tsv | awk -F "\t" '{if($4>0 && $5>0 && $7>0 && $8>0 && $11>0)print $0}' | awk -F "\t" '{if($17==0 && $18==0)print $0}' |  awk -F "\t" '{if($20==0 && $13==0 && $15==0 && $14==0 && $3==0 && $6==0 && $21==0 ) print $0}' | awk -F "\t" '{if($22==0 && $23==0 && $12==0 ) print $0}' | awk -F "\t" '{if($2==0 && $19==0 && $16==0 ) print $0}' | cut -f 1 > 04.Crassostrea

awk -F "\t" '{if($9>0 && $10>0)print $0}'  Orthogroups.GeneCount.tsv | awk -F "\t" '{if($4==0 && $5==0 && $7==0 && $8==0 && $11==0)print $0}' | awk -F "\t" '{if($17==0 && $18==0)print $0}' |  awk -F "\t" '{if($20==0 && $13==0 && $15==0 && $14==0 && $3==0 && $6==0 && $21==0 ) print $0}' | awk -F "\t" '{if($22==0 && $23==0 && $12==0 ) print $0}' | awk -F "\t" '{if($2==0 && $19==0 && $16==0 ) print $0}' | cut -f 1 >  05.Csikamea

awk -F "\t" '{if($9==0 && $10>0)print $0}'  Orthogroups.GeneCount.tsv | awk -F "\t" '{if($4==0 && $5==0 && $7==0 && $8==0 && $11==0)print $0}' | awk -F "\t" '{if($17==0 && $18==0)print $0}' |  awk -F "\t" '{if($20==0 && $13==0 && $15==0 && $14==0 && $3==0 && $6==0 && $21==0 ) print $0}' | awk -F "\t" '{if($22==0 && $23==0 && $12==0 ) print $0}' | awk -F "\t" '{if($2==0 && $19==0 && $16==0 ) print $0}' | cut -f 1 >  06.CsikameaB

awk -F "\t" '{if($9>0 && $10==0)print $0}'  Orthogroups.GeneCount.tsv | awk -F "\t" '{if($4==0 && $5==0 && $7==0 && $8==0 && $11==0)print $0}' | awk -F "\t" '{if($17==0 && $18==0)print $0}' |  awk -F "\t" '{if($20==0 && $13==0 && $15==0 && $14==0 && $3==0 && $6==0 && $21==0 ) print $0}' | awk -F "\t" '{if($22==0 && $23==0 && $12==0 ) print $0}' | awk -F "\t" '{if($2==0 && $19==0 && $16==0 ) print $0}' | cut -f 1 >  06.CsikameaA

grep -Ff 06.CsikameaB -w Orthogroups.GeneCount.tsv | \
awk -F "\t" '{
    for (i = 1; i <= 23; i++) {
        sums[i] += $i
    }
} END {
    for (i = 1; i <= 23; i++) {
        print sums[i]
    }
}'

awk -F "\t" '{print $9}' Orthogroups_UnassignedGenes.tsv | grep CsiA | sed 's/, /\t/g' | wc -l
grep -Ff 01.Bivalvia -w Orthogroups.tsv  | awk -F "\t" '{print $9}'| sed 's/, /\n/g' | awk '{print "Bivalvia""\t"$0}' > tmp.Bivalvia.csiA
grep -Ff 02.Pteriomorphia -w Orthogroups.tsv  | awk -F "\t" '{print $9}'| sed 's/, /\n/g' | awk '{print "Pteriomorphia""\t"$0}' > tmp.Pteriomorphia.csiA
grep -Ff 03.Osteridae -w Orthogroups.tsv  | awk -F "\t" '{print $9}'| sed 's/, /\n/g' | awk '{print "Osteridae""\t"$0}' > tmp.Osteridae.csiA
grep -Ff 04.Crassostrea -w Orthogroups.tsv  | awk -F "\t" '{print $9}'| sed 's/, /\n/g' | awk '{print "Crassostrea""\t"$0}' > tmp.Crassostrea.csiA
grep -Ff 05.Csikamea -w Orthogroups.tsv  | awk -F "\t" '{print $9}'| sed 's/, /\n/g' | awk '{print "Csikamea""\t"$0}' > tmp.Csikamea.csiA
grep -Ff 06.CsikameaA -w Orthogroups.tsv  | awk -F "\t" '{print $9}'| sed 's/, /\n/g' | awk '{print "CsikameaA""\t"$0}' > tmp1.CsikameaA.csiA
awk -F "\t" '{print $9}' Orthogroups_UnassignedGenes.tsv | grep CsiA | sed 's/, /\t/g' | awk '{print "CsikameaA""\t"$0}' > tmp2.CsikameaA.csiA

grep -Ff 01.Bivalvia -w Orthogroups.tsv  | awk -F "\t" '{print $10}'| sed 's/, /\n/g' | awk '{print "Bivalvia""\t"$0}' > tmp.Bivalvia.csiB
grep -Ff 02.Pteriomorphia -w Orthogroups.tsv  | awk -F "\t" '{print $10}'| sed 's/, /\n/g' | awk '{print "Pteriomorphia""\t"$0}' > tmp.Pteriomorphia.csiB
grep -Ff 03.Osteridae -w Orthogroups.tsv  | awk -F "\t" '{print $10}'| sed 's/, /\n/g' | awk '{print "Osteridae""\t"$0}' > tmp.Osteridae.csiB
grep -Ff 04.Crassostrea -w Orthogroups.tsv  | awk -F "\t" '{print $10}'| sed 's/, /\n/g' | awk '{print "Crassostrea""\t"$0}' > tmp.Crassostrea.csiB
grep -Ff 05.Csikamea -w Orthogroups.tsv  | awk -F "\t" '{print $10}'| sed 's/, /\n/g' | awk '{print "Csikamea""\t"$0}' > tmp.Csikamea.csiB
grep -Ff 06.CsikameaB -w Orthogroups.tsv  | awk -F "\t" '{print $10}'| sed 's/, /\n/g' | awk '{print "CsikameaA""\t"$0}' > tmp1.CsikameaB.csiB
awk -F "\t" '{print $10}' Orthogroups_UnassignedGenes.tsv | grep CsiB | sed 's/, /\t/g' | awk '{print "CsikameaA""\t"$0}' > tmp2.CsikameaB.csiB


#######################################################################################################################################


Identification of allelic variations
1. WGDI+BLST+EGGNOG
2. sequence similar compare (needle)
3. expr compare (STAR+stringtie)

###############################################
## ALLELE pipline
###############################################

#### BLAST
diamond makedb --in HAPA.pep.fa --db HAPA
diamond makedb --in HAPB.pep.fa --db HAPB

diamond blastp --threads 10 --db HAPB.dmnd --query HAPA.pep.fa --outfmt 6 --evalue 1e-10 --max-target-seqs 20 --out HAPA2HAPB.outblast
diamond blastp --threads 10 --db HAPA.dmnd --query HAPB.pep.fa --outfmt 6 --evalue 1e-10 --max-target-seqs 20 --out HAPB2HAPA.outblast

cut -f 1,2,12 HAPA2HAPB.outblast > HAPA2HAPB.score
cut -f 1,2,12 HAPB2HAPA.outblast > HAPB2HAPA.score
~/tools/genetribe/genetribe RBH -a HAPA2HAPB.score -b HAPB2HAPA.score > allele2blast.final.out

######  wgdi
#### Data ready
mamba  create -c bioconda -c conda-forge -n wgdi wgdi
conda activate wgdi

python .generate_conf.py -p hp1 HAP1.fasta HAPA.final.longest.gff
python 01.generate_conf.py -p hp2 HAP2.fasta HAPB.final.longest.gff

awk '{print $7"\t"$2}' hp1.gff > chang.id.hp1
awk '{print $7"\t"$2}' hp2.gff > chang.id.hp2

while read i ; do diamond makedb --in $i'.pep.fa' --db $i --threads 10; done < list
while read i ; do diamond blastp --threads 10 --db CskiameaA.dmnd --query $i'.pep.fa' --outfmt 6 --evalue 1e-5 --max-target-seqs 40 --out 'csi2'$i'.wgdi.blastp.txt' ; done < list
while read i ; do diamond blastp --threads 10 --db $i'.dmnd' --query CskiameaA.pep.fa --outfmt 6 --evalue 1e-5 --max-target-seqs 40 --out $i'2csiA.wgdi.blastp.txt' ; done < list

#### wgdi
### first get the Dotplot && block
[dotplot]
gff1 = 00.data/hp1.gff
gff2 = 00.data/hp2.gff
lens1 = 00.data/hp1.order.len
lens2 = 00.data/hp2.len
blast = 00.data/csi.wgdi.blastp.txt
blast_reverse = false
genome1_name =  C. skiamea(A)
genome2_name =  C. skiamea(B)
multiple  = 1
score = 100
evalue = 1e-5
repeat_number = 10
position = order
ancestor_left = none
ancestor_top = none
markersize = 0.5
figsize = 10,10
savefig = csi.wgdi.dot.pdf
```
wgdi -d 01.csi.wgdi.conf
```

[collinearity]
gff1 = 00.data/hp1.gff
gff2 = 00.data/hp2.gff
lens1 = 00.data/hp1.len
lens2 = 00.data/hp2.len
blast = 00.data/csi.wgdi.blastp.txt
blast_reverse = false
multiple  = 1
process = 10
evalue = 1e-5
score = 100
grading = 50,40,25
mg = 40,40
pvalue = 0.2
repeat_number = 10
positon = order
savefile = csi.wgdi.collinearity.txt
```
wgdi -icl 01.csi.wgdi.conf
```


###### Allele sequence similarity score -- needle
#### Data ready
mamba install -c bioconda emboss seqkit
mamba install -c conda-forge libiconv

cut -f 1,2 06.allele.pair | sed 's/\t/\n/g' > 07.needle/00.allele.id
seqkit grep -f 00.allele.id pep.fa > 01.allele.pep.fa

cut -f 1,2 ../06.allele.pair | awk '{print$0"\t""gene"NR}' > 02.homo

cut -f 1,3 02.homo > tmp1
cut -f 2,3 02.homo > tmp2
cat tmp1 tmp2 > 03.replace.id

seqkit replace -p '(.+)' -r '{kv}' -k 03.replace.id 01.allele.pep.fa > 04.allele.pep.fa.change

seqkit split --by-id --id-regexp "\[(.+)\]" 04.allele.pep.fa.change

for file in `ls *.change`
do
     newFile=`echo $file | sed 's/04.allele.pep.fa.part_/allele_/' | sed 's/\.change//'`
     mv $file $newFile
done

### compare by needle
while read i ; do needle -asequence '04.allele.pep.fa.change.split/'$i -bsequence '04.allele.pep.fa.change.split/'$i -gapopen 10 -gapextend 0.5 -outfile '05.result/'$i'.needle' ; done < list
grep "Similarity" *.needle > 05.needle_result

tac 05.needle_result | sed 's/\.needle:# Similarity://g' | awk '!i[$1]++' > 05.needle_result
cut -f 2 -d "(" 05.needle_result | cut -f 1 -d ")" | sed 's/%//g' > 06.tmp.needle

### gene location & length
awk '{print $9"\t"$1"\t"$4"\t"$5"\t"($5-$4+1)}' 00.all.gene.gff > 00.all.gene.stats
### pep lengh
seqkit fx2tab -l -n  01.allele.pep.fa > 01.allele.pep.table
### similiarity pep/cds/promoter
### Ka and Ks values were calculated by WGDI


##### gene function annotations --- EGGNOG
nohup emapper.py -i CskiameaA.pep.fa --output 01.CsiA -m diamond -d euk --data_dir /public1/db/eggnog-mapper-data/ --cpu 10 &
nohup emapper.py -i CskiameaB.pep.fa --output 01.CsiB -m diamond -d euk --data_dir /public1/db/eggnog-mapper-data/ --cpu 10 &

## grep -v "#" 01.CsiA.emapper.annotations > ../02.csiA.anno
## seqkit fx2tab -l -n CskiameaA.pep.fa | awk '{print $1"\t"$2}' > ../03.csiA.tab

##### using R
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
options(stringsAsFactors = F)

readxl::read_excel('ALLELE.xlsx') -> ALLEL
readxl::read_excel('ALLELE.xlsx') -> ALLEL

csiA_emapper <- read_delim('02.csiA.anno',"\t", escape_double = FALSE, trim_ws = TRUE,col_names = F)
csiB_emapper <- read_delim('02.csiB.anno',"\t", escape_double = FALSE, trim_ws = TRUE,col_names = F)

csiA_emapper <- csiA_emapper[,c(1,2,9,6)]
csiB_emapper <- csiB_emapper[,c(1,2,9,6)]

colnames(csiA_emapper) <- c("CsiA","ortholog_cgi","name_cgi","lvl_cgi")
colnames(csiB_emapper) <- c("CsiB","ortholog_csi","name_csi","lvl_csi")

merge(ALLEL,csiA_emapper,by="CsiA", all.x=T) -> ALLEL
merge(ALLEL,csiB_emapper,by="CsiB", all.x=T) -> ALLEL

#### filtered by function annotations and proteion length
data1 <- ALLEL
data1[is.na(data1)] = "-"

data2 <- data1[(data1$ortholog_csi==data1$ortholog_cgi) &  data1$ortholog_csi!="-",]
data3 <- data1[(data1$name_csi == data1$name_cgi &  data1$name_cgi!="-" ),]
data4 <- data1[data1$Apep_len>data1$Bpep_len*0.75 & data1$Apep_len*0.75 < data1$Bpep_len & (data1$ortholog_csi=="-" & data1$ortholog_cgi=="-" ),]
data5 <- rbind(data2,data3,data4)
data5<-data5[!duplicated(data5),]

data5 <- data5[data5$lvl_cgi=="33208|Metazoa" | data5$lvl_cgi=="2759|Eukaryota" | data5$lvl_cgi=="-", ]

write.csv(data5,'ALLELE.final.csv')
save.image('data.Rdata')

#### Finally, the alleles with protein sequences exhibiting more than 50% divergence were filtered. (the result of seq divergence were calculated by needle)
#### the Finally result can see in the Allele.xlsx
##################################################################################################################

allele expression
1. STAR + stringtie

### 00.trimmed and filtered
while read i ; do fastp -i $i'_1.fq.gz' -I $i'_2.fq.gz' -o $i'cleaned_1.fq.gz' -O $i'cleaned_2.fq.gz' -l 36 -q 20 -w 12 ; done < list

### 01.reads_mapping
STAR v2.7.9a was used with default parameters for the read mapping on the C.gigas reference genome with GTF file, by the Multi-sample 2-pass method.
1. Run 1st mapping pass for all samples with ”usual” parameters. Using annotations is recommended either a the genome generation step, or mapping step.
2. Run 2nd mapping pass for all samples , listing SJ.out.tab files from all samples in
--sjdbFileChrStartEnd /path/to/sj1.tab /path/to/sj2.tab ....

# 01./1/ star 1-pass index
cat HAP1.fasta HAP2.fasta > genome.fa
cat /public1/node3_liushk/tmp/03.braker_result/CskiameaA.longest.gff  /public1/node3_liushk/tmp/03.braker_result/CskiameaB.longest.gff > Cskiamea.d.gff
agat_convert_sp_gff2gtf.pl --gff Cskiamea.d.gff -o Cskiamea.d.gtf #genomic
STAR --runMode genomeGenerate --runThreadN 16 --genomeDir 00.star_index/ --genomeFastaFiles genome.fa --sjdbGTFfile Cskiamea.d.gtf --sjdbOverhang 149 #genome

# 01./2/ 1-pass mapping with indexed genome
while read i ; do STAR --genomeDir 00.star_index/ --runThreadN 20  --readFilesIn '01.fastp/'$i'_1.cleaned.fq.gz'  '01.fastp/'$i'_2.cleaned.fq.gz'  --readFilesCommand zcat --outSAMunmapped Within --outFileNamePrefix '02.1-pass/'$i  > '02.1-pass/log2map.'$i 2>&1  ; done < 01.fastp/list
cat *.tab | awk '($5 > 0 && $7 > 2 && $6==0)' | cut -f1-6 | sort | uniq > SJ.filtered.tab

# 2-pass mapping with SJ.out.tab files
while read i ; do STAR --genomeDir 00.star_index/ --runThreadN 30 --sjdbFileChrStartEnd SJ.filtered.tab --readFilesCommand zcat \
--readFilesIn '01.fastp/'$i'_1.cleaned.fq.gz' '01.fastp/'$i'_2.cleaned.fq.gz' \
--outSAMmapqUnique 255 --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 30 \
--outFileNamePrefix '03.2-pass/'$i  > '03.2-pass/log2map.'$i 2>&1 ; done <  01.fastp/list

# stringtie_expression
cat CskiameaA.longest.gff CskiameaB.longest.gff > Cskiamea.d.gff
stringtie -p 20 -G Cskiamea.d.gff  -e -B -o 04.string/TD4/transcripts.gtf -A 04.string/TD4/gene_abundances.tsv 03.2-pass/TD4Aligned.sortedByCoord.out.bam
stringtie -p 20 -G Cskiamea.d.gff  -e -B -o 04.string/TG4/transcripts.gtf -A 04.string/TG4/gene_abundances.tsv 03.2-pass/TG4Aligned.sortedByCoord.out.bam
stringtie -p 20 -G Cskiamea.d.gff  -e -B -o 04.string/TL4/transcripts.gtf -A 04.string/TL4/gene_abundances.tsv 03.2-pass/TL4Aligned.sortedByCoord.out.bam
stringtie_expression_matrix.pl --expression_metric=TPM --result_dirs='TD4,TG4,TL4' --transcript_matrix_file=transcript_tpm_all_samples.tsv --gene_matrix_file=gene_tpm_all_samples.tsv



### in R ## Fig4d

library(cowplot)
library(dplyr)
library(readr)
library(ggplot2)
library(plyr)
library(ggthemes)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(scales)
library(gghalves)
library(pheatmap)
library(ggforce)

read.csv('Allele.csv') -> aa
read.table('02.genelist.hds.csi') -> hapA
read.table('02.genelist.hds.csi.hapB') -> hapB
colnames(hapA) <- "CsiA.gene.ID"
colnames(hapB) <- "CsiB.gene.ID"
merge(hapA,aa,by="CsiA.gene.ID") -> bb
merge(hapB,bb,by="CsiB.gene.ID") -> bb



aa[aa$ks_NG86>0,] -> xx
mean(xx$ka_NG86/xx$ks_NG86)
mean(aa$pep_Similar)

read.table('02.genelist.sv.csi') -> hapA
read.table('02.genelist.sv.csi.hapB') -> hapB
colnames(hapA) <- "CsiA.gene.ID"
colnames(hapB) <- "CsiB.gene.ID"
merge(hapA,aa,by="CsiA.gene.ID") -> cc
merge(hapB,cc,by="CsiB.gene.ID") -> cc


read.table('02.genelist.synal.csi') -> hapA
read.table('02.genelist.synal.csi.hapB') -> hapB
colnames(hapA) <- "CsiA.gene.ID"
colnames(hapB) <- "CsiB.gene.ID"
merge(hapA,aa,by="CsiA.gene.ID") -> dd
merge(hapB,dd,by="CsiB.gene.ID") -> dd

## gene in hapA is 30052 / gene in hapB is 29609
## gene in hapA hds is  11105 / gene in hapB hds is 10675
### 3832

## gene in hapA SYN is  4266 / gene in hapB SYN is 4227
### 3136

## gene in hapA SV is  570 / gene in hapB SV is 632
### 141

tmp1.df <- data.frame("Group"="HDS","value"=bb$pep_Similar)
tmp2.df <- data.frame("Group"="SV","value"=cc$pep_Similar)
tmp3.df <- data.frame("Group"="SY","value"=dd$pep_Similar)
rbind(tmp1.df,tmp2.df,tmp3.df )->df1

source("geom_flat_violin.R")

df1$Group <- factor(df1$Group,levels = c("SY","SV","HDS"))

summary_simdat <- summarySE(df1, measurevar = "value",
                            groupvars = c("Group"))
comparisons <- list(c("SY", "HDS"),c("SV", "HDS"),c("SV", "SY"))


P1 <- ggplot(df1,aes(x=Group,y=value, fill = Group, colour = Group))+
    geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2)+
    geom_point(position = position_jitter(width = .15), size = .25)+
    geom_boxplot(aes(x = Group, y = value, fill = Group),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
    ylab('Score')+xlab('Group')+ theme_cowplot() + guides(fill = FALSE, colour = FALSE) +
    scale_colour_brewer(palette = "Dark2")+
    scale_fill_brewer(palette = "Dark2") + stat_compare_means(comparisons = comparisons,label = "p.signif",
                                                              method = "wilcox.test",hide.ns = F) + theme(legend.position = "none")


pdf("4d.pdf", width=8, height=4)
P1 + facet_zoom(ylim = c(95, 100),zoom.size = 0.8) + scale_fill_manual(values = c('#0d898a','#e18283','#5494cc'))+scale_color_manual(values = c('#0d898a','#e18283','#5494cc'))
dev.off()

###### Fig4e
result <- aa
result$lgTD <- (log2(result$CsiA_TD+1)-log2(result$CsiB_TD+1))
result$lgTG <- (log2(result$CsiA_TG+1)-log2(result$CsiB_TG+1))
result$lgTL <- (log2(result$CsiA_TL+1)-log2(result$CsiB_TL+1))
result1 <- result
subset(result1,abs(result1$lgTD) >1 |abs(result1$lgTL) >1 | abs(result1$lgTG) >1 ) -> result2

bk<- c(-12,-1,1,12)
P<-pheatmap (result2[,c(32:34)],breaks = bk,cellwidth = 14,cellheight = 0.08,show_rownames = F,cluster_cols = F,color =  rev(RColorBrewer::brewer.pal(3, name = "PiYG")))
ggsave(P,file="Fig4e.pdf", width=4, height=10)


#### Fig4d

merge(result2,hapA,by="CsiA.gene.ID") -> result2_hds
merge(result2_hds,hapB,by="CsiB.gene.ID") -> result2_hds



# Load the necessary libraries
library(ggplot2)
library(reshape2)

# Create the data frame
read.table('Data3.txt',header = 1) -> aa

# Reshape the data frame to long format
data_long <- melt(aa, id.vars = "GROUP", variable.name = "Category", value.name = "Value")

# Create the bar plot
P3 <- ggplot(data_long, aes(x = GROUP, y = Value, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Group", y = "Value") + theme_minimal() + theme_cowplot()  + guides(fill = FALSE, colour = FALSE)

P4 <- P3+scale_fill_manual(values =c('#1B9E77','#E7298A'))+scale_color_manual(values = c('#1B9E77','#E7298A'))
pdf("Fig4f.pdf", width=6, height=4)
P4
dev.off()
