#### /home/shicy/0922/haplotype/
# grep NOTAL syri.out | cut -f 6-8 | grep -v "-" | awk '{sum += $3 - $2+1} END {print sum}'
# awk '{sum += $3 - $2} END {print sum}'
grep NOTAL syri.out | cut -f 1-3 | grep -v "-" | awk '{print $1"\t"$2-1"\t"$3}' > syri.not.bed
bedtools sort -i syri.not.bed | awk '{print "A"$0}' > hapA.NOTAL.bed
grep NOTAL syri.out | cut -f 6-8 | grep -v "-" | awk '{print $1"\t"$2-1"\t"$3}' > syri.not.bed
bedtools sort -i syri.not.bed | awk '{print "B"$0}' > hapB.NOTAL.bed

grep HDR syri.out | cut -f 1-3 | grep -v "-" | awk '{print $1"\t"$2-1"\t"$3}' > syri.not.bed
bedtools sort -i syri.not.bed | awk '{print "A"$0}' > hapA.hdr.bed
grep HDR syri.out | cut -f 6-8 | grep -v "-" | awk '{print $1"\t"$2-1"\t"$3}' > syri.not.bed
bedtools sort -i syri.not.bed | awk '{print "B"$0}' > hapB.hdr.bed

cat hapA.NOTAL.bed hapA.hdr.bed  | sort -k1,1 -k2,2n -k3,3n  > hapA.HDS.bed
cat hapB.NOTAL.bed hapB.hdr.bed  | sort -k1,1 -k2,2n -k3,3n  > hapB.HDS.bed
cat hapA.HDS.bed hapB.HDS.bed > HDS.bed

### Repeat state
grep -v "Simple_repeat" pfuV4.1HapA_ref_genome.fasta.out.gff | grep -v "Low_complexity"  | awk '{print "A"$1"\t"$4-1"\t"$5}' | grep -v "#" > repeat.hap1.bed
grep -v "Simple_repeat" pfuV4.1HapB_alt_genome.fasta.out.gff | grep -v "Low_complexity"  | awk '{print "B"$1"\t"$4-1"\t"$5}' | grep -v "#" > repeat.hap2.bed
cat repeat.hap1.bed repeat.hap2.bed > repeat.bed

grep -v "Simple_repeat" hap1.chr.fa.out | grep -v "Low_complexity" | grep -v "#" | awk '{print "A"$5"\t"$6-1"\t"$7}' | grep Chr >  repeat.hap1.bed
grep -v "Simple_repeat" hap2.chr.fa.out | grep -v "Low_complexity" | grep -v "#" | awk '{print "B"$5"\t"$6-1"\t"$7}' | grep Chr >  repeat.hap2.bed
cat repeat.hap1.bed repeat.hap2.bed > repeat.bed

# samtools faidx 01.genome.fa
awk '{print "A"$0}' Afixchr.ref.filtered.fa.fai > HapA.fai
awk '{print "B"$0}' Afixchr.qry.filtered.fa.fai > HapB.fai
cat HapA.fai HapB.fai | awk '{print "chr - "$1" "$1" 0 "$2" genome1"}'  > karyotype.txt
cut -d ' ' -f 3,6 karyotype.txt | tr ' ' '\t' > tmp.genome
bedtools makewindows -g tmp.genome -w 500000 > tmp.windows
bedtools coverage -a tmp.windows -b 00.syri/HDS.bed | awk '{print $1,$2,$3,$7}' > HDS.density
bedtools coverage -a tmp.windows -b 00.repeat/repeat.bed | awk '{print $1,$2,$3,$7}' > repeat.density

########################## in R
library(ggpubr)

notalign_annotation <- read.table("HDS.density",header = F,stringsAsFactors = F)
colnames(notalign_annotation)<- c("Chr","Start","End","Value")
repeat_annotation <- read.table("repeat.density",header = F,stringsAsFactors = F)
colnames(repeat_annotation)<- c("Chr","Start","End","Value")

aa <- notalign_annotation
bb <- repeat_annotation

aa$ID <- paste(aa$Chr,aa$Start,aa$End,sep = "-")
bb$ID <- paste(bb$Chr,bb$Start,bb$End,sep = "-")
tmp1.df <- data.frame("ID"=aa$ID,"notalign.value"=aa$Value)
tmp2.df <- data.frame("ID"=bb$ID,"repeat.value"=bb$Value)

merge(tmp1.df,tmp2.df,by="ID") -> cor.df

## 画图
P1 <- ggscatter(cor.df, x="notalign.value", y="repeat.value", size = 1.5,
                add="reg.line", color="#449945", conf.int=T, linewidth=5,
                cor.coef = TRUE, cor.coeff.args = list(method = "pearson")) +
    labs(x="The distribution of not-align region",
         y="The distribution of TEs") +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_classic() +
    theme(axis.ticks = element_blank(),
          axis.text.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=12, face="plain"),
          axis.title.y=element_text(colour="black", size=12, face="plain"),
          axis.title.x=element_text(colour="black", size=12, face="plain"),
          panel.border = element_rect(fill=NA, color="black", size=1.5, linetype="solid"),
          panel.grid.major = element_line(size=1))

ggsave(P1, file="notalign2TEs.pdf", width=5, height=4.8)



bedtools coverage -a 00.syri/HDS.bed -b   00.repeat/repeat.bed > 01.repeatInHDS


#########################################

library(tidyverse)
library(ggplot2)
library(ggmagnify)
library(tidyverse)
library(ggplot2)
library(ggmagnify)

read.table('01.repeatInHDS.csi') -> bb
data_summary <- bb %>%
    group_by(V6) %>%
    summarize(Mean_Score = mean(V7))
p<-ggplot(data_summary, aes(x = V6, y = Mean_Score)) +
    geom_point(color = "blue", size = 0.2) +
    labs(title = "", x = "HDS length (kb)", y = "Average proportion of TEs") +
    theme_minimal()  + geom_smooth(method = "loess", color = "red") + theme_classic() +
    theme(axis.ticks = element_blank(),
          axis.text.x=element_text(colour="black",size=15),
          axis.text.y=element_text(colour="black",size=15,face="plain"),
          axis.title.y=element_text(colour="black",size = 15,face="plain"),
          axis.title.x=element_text(colour="black",size = 15,face="plain"),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"),
          panel.grid.major = element_line(size=1)) + scale_x_continuous(labels = function(x) scales::label_number(scale = 0.001)(x))

target = c(0,80000,0.25,0.75)
insert = c(600000,3000000,0,0.6)
#p1就是大图，后面就是指放大边框颜色、粗细，范围
p1 <-p + geom_magnify(from = target,
                      to = insert,
                      proj ="facing",
                      colour = "black",
                      linewidth = 0.5,
                      axes = TRUE)
pdf(file = "tmp1.pdf", width=5, height=4.8)
p1
dev.off( )



library(tidyverse)
library(ggplot2)
library(ggmagnify)
library(tidyverse)
library(ggplot2)
library(ggmagnify)

read.table('01.repeatInHDS.csi') -> bb
data_summary <- bb %>%
    group_by(V6) %>%
    summarize(Mean_Score = mean(V7))

p<-ggplot(data_summary, aes(x = V6, y = Mean_Score)) +
    geom_point(color = "blue", size = 0.2) +
    labs(title = "", x = "HDS length (kb)", y = "Average proportion of TEs") +
    theme_minimal()  + geom_smooth(method = "loess", color = "red") + theme_classic() +
    theme(axis.ticks = element_blank(),
          axis.text.x=element_text(colour="black",size=15),
          axis.text.y=element_text(colour="black",size=15,face="plain"),
          axis.title.y=element_text(colour="black",size = 15,face="plain"),
          axis.title.x=element_text(colour="black",size = 15,face="plain"),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"),
          panel.grid.major = element_line(size=1)) + scale_x_continuous(labels = function(x) scales::label_number(scale = 0.001)(x))
target = c(0,40000,0.25,0.75)
insert = c(200000,600000,0,0.6)
#p1就是大图，后面就是指放大边框颜色、粗细，范围
p1 <-p + geom_magnify(from = target,
                      to = insert,
                      proj ="facing",
                      colour = "black",
                      linewidth = 0.5,
                      axes = TRUE)
pdf(file = "csi.pdf", width=5, height=4.8)
p1
dev.off( )


library(tidyverse)
library(ggplot2)
library(ggmagnify)
library(tidyverse)
library(ggplot2)
library(ggmagnify)

read.table('01.repeatInHDS.csi') -> bb
data_summary <- bb %>%
    group_by(V6) %>%
    summarize(Mean_Score = mean(V7))

p<-ggplot(data_summary, aes(x = V6, y = Mean_Score)) +
    geom_point(color = "blue", size = 0.2) +
    labs(title = "", x = "HDS length (kb)", y = "Average proportion of TEs") +
    theme_minimal()  + geom_smooth(method = "loess", color = "red") + theme_classic() +
    theme(axis.ticks = element_blank(),
          axis.text.x=element_text(colour="black",size=15),
          axis.text.y=element_text(colour="black",size=15,face="plain"),
          axis.title.y=element_text(colour="black",size = 15,face="plain"),
          axis.title.x=element_text(colour="black",size = 15,face="plain"),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"),
          panel.grid.major = element_line(size=1)) + scale_x_continuous(labels = function(x) scales::label_number(scale = 0.001)(x))
target = c(0,40000,0.25,0.75)
insert = c(200000,600000,0,0.6)
p1 <-p + geom_magnify(from = target,
                      to = insert,
                      proj ="facing",
                      colour = "black",
                      linewidth = 0.5,
                      axes = TRUE)
pdf(file = "csi.pdf", width=5, height=4.8)
p1
dev.off( )

#########################################

read.table('01.repeatInHDS.Ase') -> bb
data_summary <- bb %>%
    group_by(V6) %>%
    summarize(Mean_Score = mean(V7))

p<-ggplot(data_summary, aes(x = V6, y = Mean_Score)) +
    geom_point(color = "blue", size = 0.2) +
    labs(title = "", x = "HDS length (kb)", y = "Average proportion of TEs") +
    theme_minimal()  + geom_smooth(method = "loess", color = "red") + theme_classic() +
    theme(axis.ticks = element_blank(),
          axis.text.x=element_text(colour="black",size=15),
          axis.text.y=element_text(colour="black",size=15,face="plain"),
          axis.title.y=element_text(colour="black",size = 15,face="plain"),
          axis.title.x=element_text(colour="black",size = 15,face="plain"),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"),
          panel.grid.major = element_line(size=1)) + scale_x_continuous(labels = function(x) scales::label_number(scale = 0.001)(x))

target = c(0,80000,0.25,0.75)
insert = c(600000,3000000,0,0.6)
p1 <-p + geom_magnify(from = target,
                      to = insert,
                      proj ="facing",
                      colour = "black",
                      linewidth = 0.5,
                      axes = TRUE)
pdf(file = "ase.pdf", width=5, height=4.8)
p1
dev.off( )

#########################################

read.table('01.repeatInHDS.Pfu') -> bb
data_summary <- bb %>%
    group_by(V6) %>%
    summarize(Mean_Score = mean(V7))

p<-ggplot(data_summary, aes(x = V6, y = Mean_Score)) +
    geom_point(color = "blue", size = 0.2) +
    labs(title = "", x = "HDS length (kb)", y = "Average proportion of TEs") +
    theme_minimal()  + geom_smooth(method = "loess", color = "red") + theme_classic() +
    theme(axis.ticks = element_blank(),
          axis.text.x=element_text(colour="black",size=15),
          axis.text.y=element_text(colour="black",size=15,face="plain"),
          axis.title.y=element_text(colour="black",size = 15,face="plain"),
          axis.title.x=element_text(colour="black",size = 15,face="plain"),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"),
          panel.grid.major = element_line(size=1)) + scale_x_continuous(labels = function(x) scales::label_number(scale = 0.001)(x))

target = c(0,40000,0.25,0.75)
insert = c(300000,1200000,0,0.6)
p1 <-p + geom_magnify(from = target,
                      to = insert,
                      proj ="facing",
                      colour = "black",
                      linewidth = 0.5,
                      axes = TRUE)
pdf(file = "pfu.pdf", width=5, height=4.8)
p1
dev.off( )
#############################################

read.table('01.repeatInHDS.Mva') -> bb
data_summary <- bb %>%
    group_by(V6) %>%
    summarize(Mean_Score = mean(V7))

p<-ggplot(data_summary, aes(x = V6, y = Mean_Score)) +
    geom_point(color = "blue", size = 0.2) +
    labs(title = "", x = "HDS length (kb)", y = "Average proportion of TEs") +
    theme_minimal()  + geom_smooth(method = "loess", color = "red") + theme_classic() +
    theme(axis.ticks = element_blank(),
          axis.text.x=element_text(colour="black",size=15),
          axis.text.y=element_text(colour="black",size=15,face="plain"),
          axis.title.y=element_text(colour="black",size = 15,face="plain"),
          axis.title.x=element_text(colour="black",size = 15,face="plain"),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"),
          panel.grid.major = element_line(size=1)) + scale_x_continuous(labels = function(x) scales::label_number(scale = 0.001)(x))

target = c(0,60000,0.2,0.5)
insert = c(500000,1800000,0.5,1)
p1 <-p + geom_magnify(from = target,
                      to = insert,
                      proj ="facing",
                      colour = "black",
                      linewidth = 0.5,
                      axes = TRUE)
pdf(file = "mva.pdf", width=5, height=4.8)
p1
dev.off( )


############################################################
####################################################################################################################################################################################
read.table('01.repeatInHDS.csi') -> bb
p<- ggplot(bb, aes(x = V6)) +
    geom_density(fill = "blue", alpha = 0.5) + coord_cartesian(ylim = c(0, 0.0006)) +
    labs(title = "", x = "HDS length (kb)", y = "Density") +
    theme_classic() +  # 使用经典主题
    theme(
        axis.ticks = element_blank(),
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15, face = "plain"),
        axis.title.y = element_text(colour = "black", size = 15, face = "plain"),
        axis.title.x = element_text(colour = "black", size = 15, face = "plain"),
        panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid"),  # 保留边框
        panel.grid = element_blank(),  # 移除网格线
        plot.background = element_blank(),  # 移除背景色
        panel.background = element_blank()  # 移除面板背景色
    ) + scale_x_continuous(labels = function(x) scales::label_number(scale = 0.001)(x)) + scale_y_continuous(labels = scales::label_number(accuracy = 0.0001))

pdf(file = "average_csi.pdf", width=5, height=4.8)
p
dev.off( )

############################################################
read.table('01.repeatInHDS.Mva') -> bb
p<- ggplot(bb, aes(x = V6)) +
    geom_density(fill = "blue", alpha = 0.5) + coord_cartesian(ylim = c(0, 0.0006)) +
    labs(title = "", x = "HDS length (kb)", y = "Density") +
    theme_classic() +  # 使用经典主题
    theme(
        axis.ticks = element_blank(),
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15, face = "plain"),
        axis.title.y = element_text(colour = "black", size = 15, face = "plain"),
        axis.title.x = element_text(colour = "black", size = 15, face = "plain"),
        panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid"),  # 保留边框
        panel.grid = element_blank(),  # 移除网格线
        plot.background = element_blank(),  # 移除背景色
        panel.background = element_blank()  # 移除面板背景色
    ) + scale_x_continuous(labels = function(x) scales::label_number(scale = 0.001)(x)) + scale_y_continuous(labels = scales::label_number(accuracy = 0.0001))

pdf(file = "average_mva.pdf", width=5, height=4.8)
p
dev.off( )


#####################################################################
orthologs state
orthofinder XXXX


awk '{if($3=="gene")print "A"$0}' ../../Pfu.ref.longest.gff > Pfu.ref.gene.gff
bedtools coverage -a Csi.ref.gene.gff -b ~/0922/haplotype/csi/00.syri/hapA.HDS.bed | awk '{if($13==1) print $0}' > 01.geneInHDS.Csi


awk '{for (i=2; i<=25; i++) if ($i == 0) next} 1' Orthogroups.GeneCount.tsv  |  awk '{sum1 += $3; sum2+=$9;sum3 += $12; sum4+=$18} END {print sum1,sum2,sum3,sum4}'

awk '{if(($3+$22) < $26 && $3!=0 && $22!=0 )print $0}' Orthogroups.GeneCount.tsv |  awk '{sum1 += $3; sum2+=$22} END {print sum1,sum2}'
awk '{if(($9+$23) < $26 && $9!=0 && $23!=0 )print $0}' Orthogroups.GeneCount.tsv |  awk '{sum1 += $9; sum2+=$23} END {print sum1,sum2}'
awk '{if(($12+$24) < $26 && $12!=0 && $24!=0 )print $0}' Orthogroups.GeneCount.tsv |  awk '{sum1 += $12; sum2+=$24} END {print sum1,sum2}'
awk '{if(($18+$25) < $26 && $18!=0 && $25!=0 )print $0}' Orthogroups.GeneCount.tsv |  awk '{sum1 += $18; sum2+=$25} END {print sum1,sum2}'


awk '{for (i=2; i<=25; i++) if ($i == 0) next} 1' Orthogroups.GeneCount.tsv | cut -f 1 > list
grep -Ff list -w Orthogroups.tsv | awk -F "\t" '{print $3}'  | sed 's/, /\n/g' | cut -f 2 -d "@" | cut -f 1 -d "." | grep g > Ase.orth2.id
# grep -Ff list -w Orthogroups.tsv | awk -F "\t" '{print $9}'  | sed 's/, /\n/g' | grep G   > CsiA.orth2.id
grep -Ff list -w Orthogroups.tsv | awk -F "\t" '{print $12}'  | sed 's/, /\n/g' |   cut -f 2 -d "@" | cut -f 1 -d "." | grep g > Mva.orth2.id
grep -Ff list -w Orthogroups.tsv | awk -F "\t" '{print $18}'  | sed 's/, /\n/g' |   cut -f 2 -d "@" | cut -f 3 -d "." | grep g > Pfu.orth2.id

grep -Ff Ase.orth2.id -w 00.geneInHDS/01.geneInHDS.Ase  | wc -l




awk '{if(($3+$22) < $26 && $3!=0 && $22!=0 )print $0}' Orthogroups.GeneCount.tsv | cut -f 1 > list2Ase
grep -Ff list2Ase -w Orthogroups.tsv | awk -F "\t" '{print $3}'  | sed 's/, /\n/g' | cut -f 2 -d "@" | cut -f 1 -d "." | grep g > Ase.orth1.id

awk '{if(($9+$23) < $26 && $9!=0 && $23!=0 )print $0}' Orthogroups.GeneCount.tsv | cut -f 1 > list2Csi
grep -Ff list2Csi -w Orthogroups.tsv | awk -F "\t" '{print $9}'  | sed 's/, /\n/g' | sed 's/@/_/g' | grep G > Csi.orth1.id

awk '{if(($12+$24) < $26 && $12!=0 && $24!=0 )print $0}' Orthogroups.GeneCount.tsv | cut -f 1 > list2Mva
grep -Ff list2Mva -w Orthogroups.tsv | awk -F "\t" '{print $12}'  | sed 's/, /\n/g' | cut -f 2 -d "@" | cut -f 1 -d "." | grep g > Mva.orth1.id

awk '{if(($18+$25) < $26 && $18!=0 && $25!=0 )print $0}' Orthogroups.GeneCount.tsv | cut -f 1 > list2Pfu
grep -Ff list2Pfu -w Orthogroups.tsv | awk -F "\t" '{print $18}'  | sed 's/, /\n/g' |   cut -f 2 -d "@" | cut -f 3 -d "." | grep g > Pfu.orth1.id


awk -F "\t" '{if($9>0 && $23>0)print $0}'  Orthogroups.GeneCount.tsv | \
awk -F "\t" '{if($3==0 && $22==0)print $0}'  | \
awk -F "\t" '{if($12==0 && $24==0)print $0}' | \
awk -F "\t" '{if($18==0 && $25==0)print $0}' | \
awk -F "\t" '{if($2==0 && $4==0 && $5==0 && $6==0 && $7==0 && $8==0)print $0}' |  \
awk -F "\t" '{if($10==0 && $11==0) print $0}' | awk -F "\t" '{if($13==0 && $14==0 && $15==0 && $16==0 && $17==0 ) print $0}' | \
awk -F "\t" '{if($19==0 && $20==0 && $21==0 ) print $0}' |  awk '{sum1 += $9; sum2+=$23} END {print sum1,sum2}'


#########################################

grep -Ff <(grep -Ff Ase.orth2.id -w 00.geneInHDS/01.geneInHDS.Ase | cut -f 9| cut -f 1 -d ";" | sed 's/ID=//g' ) -w ../Ase.ref.longest.gff | awk '{if($3=="transcript")print $9}' | cut -f 9 | cut -f 1 -d ";"  | sed 's/ID=//g' > genelist2Ase



library(org.Asenhousia.eg.db)
library(purrr)
library(tidyverse)
library(clusterProfiler)
library(dplyr)
library(stringr)
options(stringsAsFactors = F)
load("kegg_info.RData")
pathway2gene <- AnnotationDbi::select(org.Asenhousia.eg.db,keys = keys(org.Asenhousia.eg.db),
columns = c("Pathway","KO")) %>% na.omit() %>% dplyr::select(Pathway, GID)
read.table('orth2hds.ase.genelist') -> genelist
ekp <- enricher(genelist$V1,
TERM2GENE = pathway2gene,
TERM2NAME = pathway2name,
pvalueCutoff = 0.05,
qvalueCutoff = 0.05,
pAdjustMethod = "BH")
res_ekp <- as.data.frame(ekp)
write.csv(res_ekp, file = "ekp_result2ase.csv",row.names = F)



#############################################################

library(org.Mvaria.eg.db)
library(purrr)
library(tidyverse)
library(clusterProfiler)
library(dplyr)
library(stringr)
options(stringsAsFactors = F)
load("kegg_info.RData")
pathway2gene <- AnnotationDbi::select(org.Mvaria.eg.db,keys = keys(org.Mvaria.eg.db),
                                      columns = c("Pathway","KO")) %>% na.omit() %>% dplyr::select(Pathway, GID)
read.table('tmp1') -> genelist
ekp <- enricher(genelist$V1,
                TERM2GENE = pathway2gene,
                TERM2NAME = pathway2name,
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                pAdjustMethod = "BH")
res_ekp <- as.data.frame(ekp)
write.csv(res_ekp, file = "00.ekp_result2Mvaria.csv",row.names = F)

#############################################################
library(org.Pfucata.eg.db)
library(purrr)
library(tidyverse)
library(clusterProfiler)
library(dplyr)
library(stringr)
options(stringsAsFactors = F)
load("kegg_info.RData")
pathway2gene <- AnnotationDbi::select(org.Pfucata.eg.db,keys = keys(org.Pfucata.eg.db),
                                      columns = c("Pathway","KO")) %>% na.omit() %>% dplyr::select(Pathway, GID)
read.table('tmp3') -> genelist
ekp <- enricher(genelist$V1,
                TERM2GENE = pathway2gene,
                TERM2NAME = pathway2name,
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                pAdjustMethod = "BH")
res_ekp <- as.data.frame(ekp)
write.csv(res_ekp, file = "00.ekp_result2Pfucata.csv",row.names = F)

#############################################################
library(org.Csikamea.eg.db)
library(purrr)
library(tidyverse)
library(clusterProfiler)
library(dplyr)
library(stringr)
options(stringsAsFactors = F)
load("kegg_info.RData")
pathway2gene <- AnnotationDbi::select(org.Csikamea.eg.db,keys = keys(org.Csikamea.eg.db),
                                      columns = c("Pathway","KO")) %>% na.omit() %>% dplyr::select(Pathway, GID)
read.table('tmp2') -> genelist
ekp <- enricher(genelist$V1,
                TERM2GENE = pathway2gene,
                TERM2NAME = pathway2name,
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                pAdjustMethod = "BH")
res_ekp <- as.data.frame(ekp)
write.csv(res_ekp, file = "00.ekp_result2Csikamea.csv",row.names = F)

#############################################################
library(org.Asenhousia.eg.db)
library(purrr)
library(tidyverse)
library(clusterProfiler)
library(dplyr)
library(stringr)
options(stringsAsFactors = F)
load("kegg_info.RData")
pathway2gene <- AnnotationDbi::select(org.Asenhousia.eg.db,keys = keys(org.Asenhousia.eg.db),
                                      columns = c("Pathway","KO")) %>% na.omit() %>% dplyr::select(Pathway, GID)
read.table('tmp') -> genelist
ekp <- enricher(genelist$V1,
                TERM2GENE = pathway2gene,
                TERM2NAME = pathway2name,
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                pAdjustMethod = "BH")
res_ekp <- as.data.frame(ekp)
write.csv(res_ekp, file = "00.ekp_result2Asenhousia.csv",row.names = F)


###################################################################

library(ggplot2)
library(forcats)
plot2kegg <- read.csv('plot/ekp_result2plot.csv')
A<- plot2kegg
A$Description <- as.factor(A$Description)
A$Description <- fct_inorder(A$Description)
A$logQ <- -log10(A$qvalue)

ggplot(A, aes(Group, Description)) +
geom_point(aes(color=logQ, size=Rich.factor))+theme_bw()+
theme(panel.grid = element_blank(),
axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
scale_color_gradient(low='#6699CC',high='#CC3333')+
labs(x=NULL,y=NULL)+guides(size=guide_legend(order=1)) -> p

p1 <- p + guides(colour = guide_colorbar(order = 1),size = guide_legend(order = 2))
p_kegg <- p1 + theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=14,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))

pdf("plot2kegg.pdf", width=10, height=8)
p_kegg
dev.off()
