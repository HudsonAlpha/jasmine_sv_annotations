Aspects of this repository were used in the analysis for Hiatt 2024 manuscript as follows:
See: https://www.medrxiv.org/content/10.1101/2024.03.22.24304633v1.full

```
# the consensus set creation from this repository was not used for the Hiatt 2024 analysis

#jasmine_probands.sh was used to create a merged set across probands
#jasmine_annotations.sh was used to combine the merged set with various population frequency annotation VCFs
#extracted only variants called in the probands set using bcftools 

# used transfer_annotations.py to transfer the annotation values of interest (eg AF, AN, AC, etc) from the population VCFs to the proband vcf
python ../scripts/transfer_annotations.py --input annotated_calledonly.vcf --annotation_json ~/jasmine_consensus/scripts/annotations.json --annotation_db ~/jasmine_consensus/pipeline/annotation.db/data.mdb --output frequency_annotated.vcf

# cleanup - set missing genotypes to reference for proper allele counting, fix unusual BND calls that cause bcftools errors
bcftools +setGT frequency_annotated.vcf -- -t . -n 0 > missing_to_ref.vcf

python ~/jasmine_consensus/scripts/clean_bnds.py -i tagged.vcf -o clean_tagged_annotated_probands.vcf.gz

# bcftools +fill tags to calculate AN/AC etc for proband set

bcftools sort -o clean_tagged_annotated_probands.vcf.gz cleaned.vcf

tabix clean_tagged_annotated_probands.vcf.gz

#Note: bcftools view -I -s <id> is important, -I keeps it from recaculating AN and AF

while read line; do sbatch -o logs/single_${line}.log --wrap "bcftools view -I -s ${line} clean_tagged_annotated_probands.vcf.gz | bcftools filter -i 'GT=\"alt\"' -o singles/${line}.vcf.gz"; done < ../pb_halbs.tx

# extract VCFs per proband

# use all_excl.sh  all_rare.sh  all_unfilt.sh  excl_ex.sh  exclusive.sh  rare_ex.sh  rare.sh unfiltered_2.sh  unfiltered.sh to perform various filtering and counting of SVs

```
