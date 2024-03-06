def extract_values_from_bed(filename):
    # Dictionary to store start and end values for each chromosome
    chromosome_values = {}

    # Read the BED file
    with open(filename, 'r') as file:
        for line in file:
            # Split the line into columns
            columns = line.strip().split('\t')

            # Extract chromosome, start, and end values
            chromosome = columns[0]
            start = int(columns[1])
            end = int(columns[2])

            # Update values for the current chromosome
            if chromosome not in chromosome_values:
                chromosome_values[chromosome] = {'start': start, 'end': end}
            else:
                chromosome_values[chromosome]['end'] = end

    # Extract the required values for each chromosome
    for chromosome, values in chromosome_values.items():
        print(f"{chromosome}\t{values['start']}\t{values['end']}")

# Replace 'breast_and_lynch_genes_exons_processed.bed' with your actual filename
extract_values_from_bed('breast_and_lynch_genes_exons_processed.bed')


