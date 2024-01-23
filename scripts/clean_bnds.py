import argparse
import cyvcf2

def clean_bnds(input_file, output_file):
    # Open the input VCF file
    vcf = cyvcf2.VCF(input_file)

    # Create a new VCF writer for the output file
    writer = cyvcf2.Writer(output_file, vcf)

    # Iterate over each variant in the input VCF
    for variant in vcf:
        # Check if the SVTYPE is "BND"
        if variant.INFO.get('SVTYPE') == 'BND':
            # set equal to the POS because removing is awkward and this is consistant with other variant types
            variant.INFO['END'] = variant.POS

        # Write the modified variant to the output VCF
        writer.write_record(variant)

    # Close the VCF writer
    writer.close()

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Clean BND variants in a VCF file')
    parser.add_argument('-i', '--input', required=True, help='Input VCF file')
    parser.add_argument('-o', '--output', required=True, help='Output VCF file')
    args = parser.parse_args()

    # Call the clean_bnds function with the specified input and output files
    clean_bnds(args.input, args.output)
