from lmdbm import Lmdb
from cyvcf2 import VCF
import os, shutil
import pprint
import json
from itertools import islice
import argparse

def tuple_type(s):
    try:
        return tuple(s.split(','))
    except:
        raise argparse.ArgumentTypeError("Tuples must be x,y")



# Create the parser
parser = argparse.ArgumentParser(description='Build an annotation table with keys from the VCF ID field and values from annotations json file..')

# Add the arguments
parser.add_argument('-o', '--output', type=str, required=True, help='The output filename')
parser.add_argument('--annotations', required=False, help='Annotations json file.')
parser.add_argument('--vcftuples', nargs='*', type=tuple_type, help='VCF files to use if not specified in annotations. Format is vcf,description.')


# Parse the arguments
args = parser.parse_args()

if args.annotations is None and args.vcftuples is None:
    raise ValueError('Must provide either --annotations or --vcftuples')

# Now you can use args.output and args.annotations in your script
if args.annotations:
    annotations = json.load(open(args.annotations))

def batched(iterable, n):
    # batched('ABCDEFG', 3) --> ABC DEF G
    # Code from rough equivalent of itertools.batched https://docs.python.org/3/library/itertools.html#itertools.batched
    if n < 1:
        raise ValueError('n must be at least one')
    it = iter(iterable)
    while batch := tuple(islice(it, n)):
        yield batch

# snakemake object is inherited
dbfile = os.path.dirname(args.output)
print(dbfile)

if args.vcftuples:
    annotations = []
    for vcf in args.vcftuples:
        annotations.append({'vcf': vcf[0], 'description': vcf[1], 'annotations': '*'})



with Lmdb.open(file=dbfile, flag='n', map_size=1024^3) as db:
    for annotationsource in annotations:
        print(annotationsource)
        vcf = VCF(annotationsource['vcf'])
        # add an item to the db that is {description}_{annotation} = serialized vcf header info for that annotation so we can look this up to modify header of output vcf later
        if annotationsource['annotations'] == '*':
            annotationsource['annotations'] = [x['ID'] for x in vcf.header_iter() if x['HeaderType'] == 'INFO']
            print(annotationsource['annotations'])
        for field in annotationsource['annotations']:
            db[f'{annotationsource["description"]}_{field}'] = json.dumps(vcf.get_header_type(field))
        for vbatch in batched(vcf, 25_000):
            #db[v.ID] = json.dumps({f'{annotationsource["description"]}_{field}': v.INFO.get(field,".") for field in annotationsource['annotations']})
            db.update(
                {
                    v.ID: json.dumps(
                        {f'{annotationsource["description"]}_{field}': v.INFO.get(field,'.') for field in annotationsource['annotations']}
                        ) for v in vbatch
                }
            )
    #     #data = {v.ID: ";".join((f'{annotation["description"]}_{k}={v.INFO.get(k,".")}' for k in annotation['annotations'])) for v in vcf}
    #     #print(data)