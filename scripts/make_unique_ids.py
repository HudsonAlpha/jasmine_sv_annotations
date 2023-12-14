from cyvcf2 import VCF, Writer
import argparse
import sys
import os
import uuid

if __name__ == '__main__':
    p = argparse.ArgumentParser(description='Make the ID column of a VCF unique and descsriptive.')
    p.add_argument('-i', '--input', dest='input_file', help='Input VCF file', required=True)
    p.add_argument('-o', '--output', dest='output_file', help='Output VCF file', required=True)
    p.add_argument('-d', '--description', dest='description', help='Description of the VCF file to put into ID column', required=True)
    p.add_argument('--keep_id', dest='keep_ids', action='store_true', help='Keep the original IDs in the INFO column appending unique info.', default=False)
    args = p.parse_args()

    inputvcf = VCF(args.input_file)
    outputvcf = Writer(args.output_file, inputvcf)
    for variant in inputvcf:
        if args.keep_ids:
            variant.ID = f"{args.description}_{variant.ID}_{uuid.uuid4().hex}"
        else:
            variant.ID = f"{args.description}_{variant.CHROM}_{variant.POS}_{variant.INFO.get('END', variant.POS)}_{variant.INFO.get('SVTYPE')}_{uuid.uuid4().hex}"
        outputvcf.write_record(variant)
    outputvcf.close()
    inputvcf.close()