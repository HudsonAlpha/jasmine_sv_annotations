#!/bin/bash
for f in singles/*.vcf.gz; do echo $(basename $f .vcf.gz); bcftools view -R /cluster/lab/gcooper/hg38/CNV_graphing/cnv_annotation_resources/refseq_exons_plus50bp_2023-03-29.bed.gz -H $f | wc -l >  singles/$(basename $f .vcf.gz).all.unfiltered.exon.txt;  done
