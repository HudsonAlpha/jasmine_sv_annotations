#! /bin/bash
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -o logs/proband_jasmine.log
jasmine --centroid_merging --min_overlap=0.65 --min_sequence_id=0.75 --output_genotypes genome_file=/cluster/projects/pacbio_gsc_gcooper/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta out_file=tmp_jasmine.vcf threads=16 outdir=tmp_jasmine file_list=proband_sv_vcfs.txt
