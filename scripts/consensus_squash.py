from cyvcf2 import VCF, Writer
from lmdbm import Lmdb
import os
import json
from time import sleep
from lmdb import InvalidParameterError

# inherit snakemake object
# Note: jasmine counts samples starting with 0


def combine_genotypes(cyvcf_genotype):
    assert len(cyvcf_genotype) == 2 # require two sample genotypes in
    for g in cyvcf_genotype:
        assert (g[2] is True) or ([g[0],g[1]] != [1,0]) # ensure phased 1/0 are marked
        for allele in g[:2]:
            assert allele <= 1 # ensure no multi-allelic calls
        if g[2] is True: # remove phasing in case of 1/0 
            if g[0] == 1 and g[1] == 0:
                g[0] = 0
                g[1] = 1
    second_allele = max(cyvcf_genotype[0][1], cyvcf_genotype[1][1]) # essentially, return a 1 if we see it
    first_allele = min(cyvcf_genotype[0][0], cyvcf_genotype[1][0]) # return a 0 if anyone has called 0/1
    if first_allele < 0:
        first_allele = 0
    final_gt = [[first_allele,second_allele,False],[first_allele,second_allele,False]]
    return final_gt


# simple combination: output_genotype = [max(cyvcf_genotype[0][0], cyvcf_genotype[1][0]), max(cyvcf_genotype[0][1], cyvcf_genotype[1][1]), False]

caller_key = snakemake.params.ordered_callers
sourcefiles = snakemake.input['source_vcfs']
dbfiles = snakemake.input['annotations_dbs'] # assume in same order defined by ordered callers
input_vcf = VCF(snakemake.input['vcf'])

attempts=0
for caller in caller_key:
    caller['source'] = sourcefiles[caller['order']]
    while attempts < 10:
        try:
            caller['dbfile'] = Lmdb.open(dbfiles[caller['order']], 'r')
            break
        except InvalidParameterError:
            print(f'Waiting for {dbfiles[caller["order"]]} to unlock... attempt {attempts}')
            attempts += 1
            sleep(10)
    check_vcf = VCF(caller['source']['vcf'])
    caller['annotations'] = [x['ID'] for x in check_vcf.header_iter() if x['HeaderType'] == 'INFO']
    for annotation in caller['annotations']:
        annokey = f'{caller["name"]}_{annotation}'
        annoinfo = json.loads(caller['dbfile'][annokey])
        del annoinfo['IDX'] # remove the IDX field that cyvcf2 tracks, won't be relevant in new vcf
        annoinfo['ID'] = annokey # rename annotation
        annoinfo['Description'] = annoinfo['Description'].replace('"', '') # clean up quotation mark hell
        annoinfo['Description'] = f'{caller["name"]} field {annoinfo["Description"]}' # add annotation source to description
        print(f'Adding to header: {annokey} {annoinfo}')
        input_vcf.add_info_to_header(annoinfo)


output_vcf = Writer(output_vcf_path, input_vcf)
for variant in input_vcf:
        variant.genotype = combine_genotypes(variant.genotypes)

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
