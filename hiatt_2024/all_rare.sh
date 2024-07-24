#!/bin/bash
for f in singles/*.vcf.gz; do echo $(basename $f .vcf.gz); bcftools view -R /cluster/lab/gcooper/hg38/CNV_graphing/cnv_annotation_resources/refseq_exons_plus50bp_2023-03-29.bed.gz -e 'AC>3 | inhouse_pbsv_AC >3 | gnomad4_POPMAX_AF >= 0.01 | hgsvc2_AF >= 0.01 | hprc_giab_pbsv_AF >= 0.01' -H $f | wc -l >  singles/$(basename $f .vcf.gz).all.rare.exon.txt;  done
