Fig 1C

## orthofinder
##

cp -r ../13.orthofinder/02.orth2tree/OrthoFinder/Results_Oct27/Single_Copy_Orthologue_Sequences/ 00.Single_Copy_Orthologue_Sequences
for i in *.fa ; do grep ">" $i | sed 's/>//g' > $i.1  ; done
ls *fa | cut -f 1 -d "." > 00.list

while read i; do seqkit grep -f 00.Single_Copy_Orthologue_Sequences/$i.fa.1 01.cds.fa > 02.Single_Copy_cds/$i.cds ; done < 00.Single_Copy_Orthologue_Sequences/00.list

split -l 80 02.Single_Copy_cds/list -a 2 -d  omm_

while read i; do singularity exec -B /public1:/public1 /public1/db/sif/omm_macse_v11.05b.sif /OMM_MACSE/S_OMM_MACSE_V11.05.sh --in_seq_file '01.csi_gene/'$i'.cds' --out_dir '03.omm_csi/00.omm/'$i --out_file_prefix $i --java_mem 8000m ; done < 03.omm_csi/omm_31

while read i; do singularity exec -B /public1:/public1 /public1/db/sif/omm_macse_v11.05b.sif /OMM_MACSE/S_OMM_MACSE_V11.05.sh --in_seq_file '03.split.gene/'$
i --out_dir '04.omm/omm_00_re/'$i --out_file_prefix $i --java_mem 8000m ; done < 04.omm/omm_00

grep -c ">" */*_final_align_NT.aln > 00.list
awk -F ":" '{if($2=="14")print$1}' 00.list > 01.list  # 6594
while read i ; do cp $i ../01.omm2paml/; done < 01.list
for i in *.aln ; do awk -F "@" '{print $1}' $i > $i'.fas' ; done
seqkit concat *.fas  >  ../single_copy.csd_msa.fasta
modeltest-ng -i single_copy.cds_msa.fasta -d nt -o cds -p 40
raxml-ng --msa single_copy.cds_msa.fasta --model GTR+I+G4 --threads 40 --bs-trees 200 --seed 1028 --all --outgroup Npompilius

trimal -in single_copy.cds_msa.fasta -out single_copy.cds_msa.phy -phylip_paml

perl extract_4d_phy.pl single_copy.cds_msa.phy single_copy.cds_msa.4d.phy
cp /public1/node3_shicy/PRJ02_Csikamea/14.tree/06.divtime/single_copy.cds_msa.4d.phy 01.single_copy.cds_msa.4d.phy
21 1
(((Acalifornica,Pcanaliculata),((Mmercenaria,(Tcrocea,Tgigas)),(((((Ode,Oed),((Cni,(Car,((Cgi,Can),CsiA))),Cvi)),Pfucata),(Mvirgata,MytCali)),((chfa,Myessoensis),(apu,Pmaximus))))),Npompilius)'@5.48';
baseml baseml.ctl


grep "+" mlb
0.215541 +- 0.002149

21 1
(((Acalifornica,Pcanaliculata)'B(4.061,4.215)',((Mmercenaria,(Tcrocea,Tgigas)),(((((Ode,Oed),((Cni,(Car,((Cgi,Can),CsiA))),Cvi)),Pfucata),(Mvirgata,MytCali)'B(3.179,3.409)'),((chfa,Myessoensis),(apu,Pmaximus))))),Npompilius)'B(4.8,5.594)';
mcmctree mcmctree.ctl

cp -r /public1/node3_shicy/PRJ02_Csikamea/02.haplotype_expression/05.expr2region .
cp ../01.hap.orth/OrthoFinder/Results_Oct27/Orthogroups/Orthogroups.GeneCount.tsv .

#####################
Fig1D
minimap2 -t 48 -ax asm5 --eqx ref.genome.fa qry.genome.fa > out.sam
~/shicy/fixchr/bin/fixchr -c out.sam -r ref.genome.fa -q qry.genome.fa -F S
plotsr --sr syri.out --genomes genome.txt -S 0.4 -o genome.pdf

#### HDS identity
grep NOTAL syri.out | cut -f 1-3 | grep -v "-" | awk '{print $1"\t"$2-1"\t"$3}' > syri.not.bed
bedtools sort -i syri.not.bed | awk '{print "A"$0}' > hapA.NOTAL.bed
grep NOTAL syri.out | cut -f 6-8 | grep -v "-" | awk '{print $1"\t"$2-1"\t"$3}' > syri.not.bed
bedtools sort -i syri.not.bed | awk '{print "B"$0}' > hapB.NOTAL.bed

grep HDR syri.out | cut -f 1-3 | grep -v "-" | awk '{print $1"\t"$2-1"\t"$3}' > syri.not.bed
bedtools sort -i syri.not.bed | awk '{print "A"$0}' > hapA.hdr.bed
grep HDR syri.out | cut -f 6-8 | grep -v "-" | awk '{print $1"\t"$2-1"\t"$3}' > syri.not.bed
bedtools sort -i syri.not.bed | awk '{print "B"$0}' > hapB.hdr.bed

cat hapA.NOTAL.bed hapA.hdr.bed | sort -k1,1 -k2,2n -k3,3n  > hapA.HDS.bed
cat hapB.NOTAL.bed hapB.hdr.bed | sort -k1,1 -k2,2n -k3,3n  > hapB.HDS.bed

############################################################
Fig 1e

grep NOTAL syri.out | cut -f 1-3 | grep -v "-" | awk '{print $1"\t"$2-1"\t"$3}' > syri.not.bed
bedtools sort -i syri.not.bed | awk '{print $0}' > hapA.NOTAL.bed
grep NOTAL syri.out | cut -f 6-8 | grep -v "-" | awk '{print $1"\t"$2-1"\t"$3}' > syri.not.bed
bedtools sort -i syri.not.bed | awk '{print $0}' > hapB.NOTAL.bed

grep HDR syri.out | cut -f 1-3 | grep -v "-" | awk '{print $1"\t"$2-1"\t"$3}' > syri.not.bed
bedtools sort -i syri.not.bed | awk '{print $0}' > hapA.hdr.bed
grep HDR syri.out | cut -f 6-8 | grep -v "-" | awk '{print $1"\t"$2-1"\t"$3}' > syri.not.bed
bedtools sort -i syri.not.bed | awk '{print $0}' > hapB.hdr.bed

cat hapA.NOTAL.bed hapA.hdr.bed | sort -k1,1 -k2,2n -k3,3n  > hapA.HDS.bed
cat hapB.NOTAL.bed hapB.hdr.bed | sort -k1,1 -k2,2n -k3,3n  > hapB.HDS.bed

cat hapA.HDS.bed hapB.HDS.bed > car.HDS.bed
bedtools coverage -a tmp.windows -b car.HDS.bed | awk '{print $1,$2,$3,$7}' > can.HDS.density


############################################################

library(ggpubr)

notalign_annotation <- read.table("Car.HDS.density",header = F,stringsAsFactors = F)
colnames(notalign_annotation)<- c("Chr","Start","End","Value")
repeat_annotation <- read.table("Car.repeat.density",header = F,stringsAsFactors = F)
colnames(repeat_annotation)<- c("Chr","Start","End","Value")

aa <- notalign_annotation
bb <- repeat_annotation

aa$ID <- paste(aa$Chr,aa$Start,aa$End,sep = "-")
bb$ID <- paste(bb$Chr,bb$Start,bb$End,sep = "-")

tmp1.df <- data.frame("ID"=aa$ID,"notalign.value"=aa$Value)
tmp2.df <- data.frame("ID"=bb$ID,"repeat.value"=bb$Value)

merge(tmp1.df,tmp2.df,by="ID") -> cor.df

P1 <- ggscatter(cor.df, x="notalign.value", y="repeat.value", size = 1.5,
                add="reg.line", color="#449945", conf.int=T, linewidth=5,
                cor.coef = TRUE,cor.coeff.args = list(method = "pearson")) +
    labs(  x="The distribution of not-align region",
           y="The distribution of TEs") +
    theme_classic() +
    theme(axis.ticks = element_blank(),
          axis.text.x=element_text(colour="black",size=12),
          axis.text.y=element_text(colour="black",size=12,face="plain"),
          axis.title.y=element_text(colour="black",size = 12,face="plain"),
          axis.title.x=element_text(colour="black",size = 12,face="plain"),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"),
          panel.grid.major = element_line(size=1))

ggsave(P1,file="notalign2TEs.pdf", width=5, height=4.8)



grep NOTAL syri.out | cut -f 1-3 | grep -v "-" | awk '{print $1"\t"$2-1"\t"$3}' > syri.not.bed
bedtools sort -i syri.not.bed | awk '{print "A"$0}' > hapA.NOTAL.bed
grep NOTAL syri.out | cut -f 6-8 | grep -v "-" | awk '{print $1"\t"$2-1"\t"$3}' > syri.not.bed
bedtools sort -i syri.not.bed | awk '{print "B"$0}' > hapB.NOTAL.bed

grep HDR syri.out | cut -f 1-3 | grep -v "-" | awk '{print $1"\t"$2-1"\t"$3}' > syri.not.bed
bedtools sort -i syri.not.bed | awk '{print "A"$0}' > hapA.hdr.bed
grep HDR syri.out | cut -f 6-8 | grep -v "-" | awk '{print $1"\t"$2-1"\t"$3}' > syri.not.bed
bedtools sort -i syri.not.bed | awk '{print "B"$0}' > hapB.hdr.bed

cat hapA.NOTAL.bed hapA.hdr.bed | sort -k1,1 -k2,2n -k3,3n  > hapA.HDS.bed
cat hapB.NOTAL.bed hapB.hdr.bed | sort -k1,1 -k2,2n -k3,3n  > hapB.HDS.bed


awk '{sum1 += $3-$2} END {print sum1}' hapA.hdr.bed
awk '{sum1 += $3-$2} END {print sum1}' hapB.hdr.bed
awk '{sum1 += $3-$2} END {print sum1}' hapA.NOTAL.bed
awk '{sum1 += $3-$2} END {print sum1}' hapB.NOTAL.bed
wc -l hapA.hdr.bed
wc -l hapB.hdr.bed
wc -l hapA.NOTAL.bed
wc -l hapB.NOTAL.bed
