from lmdbm import Lmdb
from cyvcf2 import VCF
import os, shutil
import pprint
import json
from itertools import islice

def batched(iterable, n):
    # batched('ABCDEFG', 3) --> ABC DEF G
    # Code from rough equivalent of itertools.batched https://docs.python.org/3/library/itertools.html#itertools.batched
    if n < 1:
        raise ValueError('n must be at least one')
    it = iter(iterable)
    while batch := tuple(islice(it, n)):
        yield batch

# snakemake object is inherited
dbfile = os.path.dirname(snakemake.output[0])
print(dbfile)

with Lmdb.open(file=dbfile, flag='n', map_size=1024^3) as db:
    for annotationsource in snakemake.params.annotations:
        print(annotationsource)
        try: # primary case, VCF is from the annotation dictionary
            vcf = VCF(annotationsource['vcf'])
        except OSError: # secondary case, the input declaration the inputs key
            vcf = VCF(snakemake.input[annotationsource['vcf']])
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