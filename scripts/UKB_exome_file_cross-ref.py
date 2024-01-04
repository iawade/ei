import argparse

def main(input_bed, output_txt, reference_file):
    # Initialize a dictionary to store BED information
    bed_d = {}

    # Read input BED file and populate bed_d dictionary
    with open(input_bed) as ipt:
        for line in ipt:
            line = line.strip().split('\t')
            chrom = int(line[0].replace('chr', '').replace('X', '23'))
            if chrom not in bed_d:
                bed_d[chrom] = []
            bed_d[chrom].append([int(line[1]), int(line[2])])

    # Sort BED coordinates for each chromosome
    for chrm in bed_d:
        bed_d[chrm] = sorted(bed_d[chrm], key=lambda x: (x[0], x[1]))

    # Open the output file for writing
    with open(output_txt, 'w') as opt:
        # Read reference file listing pVCF block coordinates
        with open(reference_file) as ipt:
            file_names = set()  # To store unique file names

            for line in ipt:
                line = line.strip().split('\t')
                if line[1] == '23':
                    file_chrom = 'X'
                else:
                    file_chrom = line[1]
                file_name = f'ukb23157_c{file_chrom}_b{line[2]}_v1.vcf.gz'
                if int(line[1]) in bed_d:
                    for exon in bed_d[int(line[1])]:
                        if int(line[3]) <= exon[0] <= int(line[4]) or \
                                int(line[3]) <= exon[1] <= int(line[4]):
                            file_names.add(file_name)

            # Write the sorted and unique file names to the output file
            for file_name in sorted(file_names):
                opt.write(f"{file_name}\n")

if __name__ == "__main__":
    # Set up argparse to read command-line arguments
    parser = argparse.ArgumentParser(description='Process BED file and generate output.')
    parser.add_argument('input_bed', help='Path to input BED file')
    parser.add_argument('output_txt', help='Path to output TXT file')
    parser.add_argument('reference_file', help='Path to reference file')
    args = parser.parse_args()

    # Call the main function with provided arguments
    main(args.input_bed, args.output_txt, args.reference_file)

