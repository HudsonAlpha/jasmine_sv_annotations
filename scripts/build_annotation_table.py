from lmdbm import Lmdb
from cyvcf2 import VCF
import os, shutil
import pprint

# snakemake object is inherited
dbfile = os.path.dirname(snakemake.output[0])
print(dbfile)
try:
    shutil.rmtree(dbfile)
except FileNotFoundError:
    pass


with Lmdb.open(dbfile, 'c') as db:
    # db['test_record'] = 'test_record'
    # print(db['test_record'])
    # print('done test record')
    for annotation in snakemake.params.annotations:
        print(annotation)
        vcf = VCF(annotation['vcf'])
        # add an item to the db that is {description}_{annotation} = serialized vcf header info for that annotation so we can look this up to modify header of output vcf later
        #print(vcf)
        #db.update({v.ID: "test" for v in vcf})
        db.update(
            {
                v.ID: ";".join(
                    (f'{annotation["description"]}_{k}={v.INFO.get(k,".")}' for k in annotation['annotations'])
                ) for v in vcf
            }
        )
    #     #data = {v.ID: ";".join((f'{annotation["description"]}_{k}={v.INFO.get(k,".")}' for k in annotation['annotations'])) for v in vcf}
    #     #print(data)