from cyvcf2 import VCF, Writer
from lmdbm import Lmdb
import os
import json

# inherit snakemake object
dbfile = os.path.dirname(snakemake.input['annotation_db'])
print(dbfile)
annotation_db = Lmdb.open(dbfile, 'r')
input_vcf_path = snakemake.input['vcf']
output_vcf_path = snakemake.output['vcf']

input_vcf = VCF(input_vcf_path)
for annotationsource in snakemake.params.annotations:
    for annotation in annotationsource['annotations']:
        annokey = f'{annotationsource["description"]}_{annotation}'
        annoinfo = json.loads(annotation_db[annokey])
        del annoinfo['IDX'] # remove the IDX field that cyvcf2 tracks, won't be relevant in new vcf
        annoinfo['ID'] = annokey # rename annotation
        annoinfo['Description'] = annoinfo['Description'].replace('"', '') # clean up quotation mark hell
        annoinfo['Description'] = f'{annotationsource["description"]} field {annoinfo["Description"]}' # add annotation source to description
        print(f'Adding to header: {annokey} {annoinfo}')
        input_vcf.add_info_to_header(annoinfo)

output_vcf = Writer(output_vcf_path, input_vcf)
for variant in input_vcf:
        for fk in variant.INFO.get('IDLIST',"").split(','):
            try:
                annos_to_add = json.loads(annotation_db[fk])
            except KeyError:
                annos_to_add = {}
            print(f'{annos_to_add} to {variant.CHROM} {variant.POS} {variant.REF} {variant.ALT}')
            for k,v in annos_to_add.items():
                variant.INFO[k] = v
        output_vcf.write_record(variant)
output_vcf.close()

