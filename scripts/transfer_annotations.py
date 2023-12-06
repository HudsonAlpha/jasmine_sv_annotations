from cyvcf2 import VCF, Writer
import argparse
import sys
import os

def filename_tuple_type(arg):
    try:
        filename, description, keys = arg.split(',', 2)
        if not os.path.isfile(filename):
            raise argparse.ArgumentTypeError(f"{filename} is not a valid file path.")
        keys = keys.split(',')
        if len(keys) == 0:
            raise argparse.ArgumentTypeError("No keys were provided.")
        return filename, description, keys
    except ValueError:
        raise argparse.ArgumentTypeError("Tuples must be in the format 'filename,description,keys'.")

def has_intersection(collection1, collection2):
    set1 = set(collection1)
    set2 = set(collection2)
    return bool(set1 & set2)  # '&' is the set intersection operator

if __name__ == "__main__":
    p = argparse.ArgumentParser(description='Transfer annotations from input VCFs to a Jasmine-merged target VCF')
    p.add_argument('-i', '--input', dest='input_file', help='Input VCF file', required=True)
    p.add_argument('-a', '--annotation', dest="annotation_tuples", metavar='filename,description,keys', type=filename_tuple_type, nargs='+', help='Tuples of filename and description and fieldlist separated by a comma. Fieldlist is a comma-separated list of annotations to transfer', required=True)
    p.add_argument('-o', '--output', dest='output_file', help='Output VCF file', required=True)
    args = p.parse_args()

    inputvcf = VCF(args.input_file)
    annotationVCFs = [list(VCF(filename)) for filename, description, keys in args.annotation_tuples]
    annotations = [keys for filename, description, keys in args.annotation_tuples]
    descriptions = [description for filename, description, keys in args.annotation_tuples]
    outputvcf = Writer(args.output_file, inputvcf) #inputvcf is used as a template (second argument)

    # build the union of all available annotation IDs.
    # Instead: Look at the support vector minus the first one (assuming first is input)
    # For each position of the support vector lookup the corresponding ID in the [position of 1 in support vector]th VCF
    # The Looking up based on the nth id in IDLIST where n is the: the current support vector position is the nth 1 in the support vector

for variant in inputvcf:
    
    suppvec = variant.INFO.get('SUPP_VEC')
    idlist = variant.INFO.get('IDLIST')
    print(variant)
    if idlist is None or suppvec is None:
        continue
    # assume first entry in suppvec represents input vcf
    originalsuppvec = suppvec
    suppvec = suppvec[1:]
    # print(range(len(suppvec)))
    # print(suppvec)
    assert(len(suppvec) == len(annotationVCFs))
    for index in range(len(suppvec)):
        entry = suppvec[index]
        if int(entry) == 1:
            rank = originalsuppvec[:index+1].count('1') # +1 because we assume first entry is input vcf
            id = idlist.split(',')[rank]
            description = descriptions[index]
            anno_list = annotations[index]
            print(f'\t\tFind {id} {anno_list} {description} {rank}')
            annotationVCF = annotationVCFs[index]
            #anno_iterator = iter(annotationVCF)
            result = [v for v in annotationVCF if v.ID == id]
            print(f'\t\t{result}')


            



    # all_annotation_ids = set([])
    # idsets = [ set([v.ID for v in annotationVCF]) for v in annotationVCFs ]
    # for idset in idsets:
    #     all_annotation_ids.update(idset)
    
    # for variant in inputvcf:
    #     component_ids = variant.INFO.get('IDLIST').split(',')
    #     if len(component_ids) == 0:
    #         continue
    #     if not has_intersection(component_ids,all_annotation_ids):
    #         continue
    #     findmatching = []
            










