import argparse
import gzip

def tally_variants(input_file, output_file):
    with gzip.open(input_file, 'rt') as infile, open(output_file, 'w') as outfile:
        # Read the header
        header_line = infile.readline().rstrip('\n')

        # Skip lines starting with "##"
        while header_line.startswith("##"):
            header_line = infile.readline().rstrip('\n')

        # Ensure the header line starts with a single "#"
        if not header_line.startswith("#"):
            raise ValueError("Invalid VCF header")

        # Split the header line into a list of column names
        columns = header_line[1:].split('\t')

        # Initialize a dictionary to store the tally for each individual
        individual_tally = {}

        # Initialize a dictionary to store the matching lines and genotypes for each individual
        individual_matches = {}

        # Process each line and update the tally for each individual
        for line in infile:
            if line.startswith("#"):
                continue  # Skip lines starting with "#"

            fields = line.rstrip('\n').split('\t')

            # Iterate through numeric columns and check for non-homref genotypes
            for eid, value in zip(columns[9:], fields[9:]):  # Start from index 9 to skip the non-numeric columns
                # Initialize the tally for the individual if not already present
                if eid not in individual_tally:
                    individual_tally[eid] = 0
                    individual_matches[eid] = {"lines": [], "genotypes": []}

                # Flag to check if any non-homref genotype is found in the line
                non_homref_found = False

                try:
                    # Extract the first, second, fourth, and fifth fields
                    matching_fields = [fields[0], fields[1], fields[3], fields[4]]

                    # Join the matching fields with "_"
                    matching_line = "_".join(matching_fields)

                    # Extract the first field separated by ":"
                    geno_value = value.split(':')[0]

                    # Check for non-homref genotype
                    if '.' not in geno_value and int(geno_value.split("/")[0]) != int(geno_value.split("/")[1]):
                        non_homref_found = True
                except (IndexError, ValueError):
                    pass  # Ignore errors if the field is not as expected

                # Update the tally and matching information if a non-homref genotype is found
                if non_homref_found:
                    individual_tally[eid] += 1
                    individual_matches[eid]["lines"].append(matching_line)
                    individual_matches[eid]["genotypes"].append(value)

        # Write the tally, matching lines, and genotypes for each individual to the output file
        outfile.write("eid\tnon_homref_count\tmatching_lines\tmatching_genotypes\n")
        for eid in individual_tally.keys():
            tally = individual_tally[eid]
            matching_lines = ";".join(individual_matches[eid]["lines"]) if tally > 1 else individual_matches[eid]["lines"][0] if individual_matches[eid]["lines"] else ""
            matching_genotypes = ";".join(individual_matches[eid]["genotypes"]) if tally > 1 else individual_matches[eid]["genotypes"][0] if individual_matches[eid]["genotypes"] else ""

            outfile.write(f'{eid}\t{tally}\t{matching_lines}\t{matching_genotypes}\n')

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Tally variants in a compressed input file.")
    parser.add_argument("input_file", help="Path to the compressed input file (e.g., input.vcf.gz)")
    parser.add_argument("output_file", help="Path to the output file with the tally for each individual")

    args = parser.parse_args()

    # Call the function with command line arguments
    tally_variants(args.input_file, args.output_file)
