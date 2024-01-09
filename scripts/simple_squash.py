from numpy import argmax
# inherit snakemake object
sample = snakemake.wildcards['sample']

with open(snakemake.input['vcf'], 'r') as input_vcf, open(snakemake.output['vcf'], 'w') as output_vcf:
    for line in input_vcf:
        line = line.strip().split('\t')
        if line[0].startswith('#CHROM'):
            output_vcf.write('\t'.join(line[:9]) + '\t' + sample + '\n')
        elif line[0].startswith('#'):
            output_vcf.write('\t'.join(line) + '\n')
        else:
            ID = line[2]
            # assume jasmine format is always GT:IS:OT:DV:DR
            # first sample is index 9
            try:
                SUPP_VEC = [x.split('=')[1] for x in line[7].split(';') if x.split('=')[0] == 'SUPP_VEC'][0]
            except IndexError:
                SUPP_VEC = '00'
            if SUPP_VEC.count('1') >= 2: #change this to "if support vector has more than one 1"
                try:
                    IDLIST = [x.split('=')[1] for x in line[7].split(';') if x.split('=')[0] == 'IDLIST'][0]
                except IndexError:
                    IDLIST = ID
                ID = IDLIST
                DV = (int(line[9].split(':')[3]),int(line[10].split(':')[3])) # Pick the genotype info with the highest DV, number supporting reads for the variant
                chosen_variant_index = argmax(DV)
            elif SUPP_VEC.count('1') == 1: # change this to if support vector only has a single 1, and determine the index of it
                chosen_variant_index = SUPP_VEC.index('1')
            else:
                chosen_variant_index = 0 # should not get here when there is a support vector but if we do, take first genotype available
                print('Warning: choising first variant index at {line[0:3]}')
            output_vcf.write('\t'.join(line[:2]) + '\t' + ID + '\t' + '\t'.join(line[3:9]) + '\t' + line[9+chosen_variant_index] + '\n')

