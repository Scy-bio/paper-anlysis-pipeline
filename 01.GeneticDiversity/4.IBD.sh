#!/bin/bash
##IBD
beagle -Xmx60g -Djava.io.tmpdir=./tmp gt=all_snp.recode.vcf.gz  out=./out_beagle nthreads=9#phasing
java -Xmx150g -jar /public/home/chcg/software/refined-ibd.17Jan20.102.jar   gt=./out_beagle.vcf.gz out=./out_beagle_phase.idb window=5e-3 trim=5e-6 lod=5 length=1e-5 

