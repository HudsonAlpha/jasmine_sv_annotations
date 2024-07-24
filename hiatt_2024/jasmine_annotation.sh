#! /bin/bash
#SBATCH -c 16
#SBATCH --mem=256G
#SBATCH -o logs/annotation_jasmine.log
jasmine --comma_filelist --centroid_merging --min_overlap=0.65 --min_sequence_id=0.75 genome_file=/cluster/projects/pacbio_gsc_gcooper/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta out_file=anno_jasmine.vcf threads=16 outdir=anno_jasmine file_list=tmp_jasmine.vcf,/cluster/projects/pacbio_gsc_gcooper/resources/sv_annotations/original/gnomad.v4.0.sv.ALL.vcf,/cluster/projects/pacbio_gsc_gcooper/resources/sv_annotations/original/hgsvc2_tagged_variants_freeze4_sv_insdel_alt.uniqueid.vcf,/cluster/projects/pacbio_gsc_gcooper/resources/sv_annotations/original/HPRC_GIAB.GRCh38.pbsv.uniqueid.vcf
