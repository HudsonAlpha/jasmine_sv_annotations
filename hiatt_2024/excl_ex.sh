#!/bin/bash
for f in singles/*.vcf.gz; do echo $(basename $f .vcf.gz); bcftools filter -R /cluster/lab/gcooper/hg38/CNV_graphing/cnv_annotation_resources/refseq_exons_plus50bp_2023-03-29.bed.gz -e 'AC>1 | inhouse_pbsv_AC>2 | gnomad4_AC>0 | hgsvc2_AC>0 | hprc_giab_pbsv_AC>0' $f | bcftools query -f '%SVTYPE\n' | sort | uniq -c > singles/$(basename $f .vcf.gz).counts.exclusive.exon.txt;  done
