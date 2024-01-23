from cyvcf2 import VCF, Writer
from lmdbm import Lmdb
import os
import json
from time import sleep
from lmdb import InvalidParameterError
import argparse
def tuple_type(s):
    try:
        return tuple(s.split(','))
    except:
        raise argparse.ArgumentTypeError("Tuples must be x,y")



parser = argparse.ArgumentParser(description='Transfer ID-keyed annotations from an annotation table to a Jasmine VCF.')
parser.add_argument('-o', '--output', type=str, required=True, help='The output filename')
parser.add_argument('-a', '--annotation_db', required=True, help='Annotation table.')
parser.add_argument('-j', '--annotation_json', required=False, help='Annotation json file.')
parser.add_argument('-i', '--input', required=True, help='Input VCF.')
parser.add_argument('--vcftuples', nargs='*', type=tuple_type, help='VCF files to use if not specified in annotations. Format is vcf,description.')
args = parser.parse_args()

if args.annotation_json is None and args.vcftuples is None:
    raise ValueError('Must provide either --annotations or --vcftuples')

if args.annotation_json:
    annotations = json.load(open(args.annotation_json))

if args.vcftuples:
    annotations = []
    for vcf in args.vcftuples:
        annotations.append({'vcf': vcf[0], 'description': vcf[1], 'annotations': '*'})

# inherit snakemake object
dbfile = os.path.dirname(args.annotation_db)
print(dbfile)
attempts = 0
while attempts < 10:
    try:
        annotation_db = Lmdb.open(dbfile, 'r')
        break
    except InvalidParameterError:
        print(f'Waiting for db to unlock... attempt {attempts}')
        attempts += 1
        sleep(10)

input_vcf_path = args.input
output_vcf_path = args.output


input_vcf = VCF(input_vcf_path)
for annotationsource in annotations:
    if annotationsource['annotations'] == '*':
        check_vcf = VCF(annotationsource['vcf'])
        annotationsource['annotations'] = [x['ID'] for x in check_vcf.header_iter() if x['HeaderType'] == 'INFO']
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
            #print(f'{annos_to_add} to {variant.CHROM} {variant.POS} {variant.REF} {variant.ALT}')
            for k,v in annos_to_add.items():
                if input_vcf.get_header_type(k)['Type'] == 'Flag' and v == '.': # This will properly remove flag annotations instead of showing FLAG='.'
                    variant.INFO[k] = False
                    continue
                try:
                    variant.INFO[k] = v
                except AttributeError: # attempt to handle strangely formatted INFO fields (tuples,etc.)
                    variant.INFO[k] = str(k)
        output_vcf.write_record(variant)
output_vcf.close()

