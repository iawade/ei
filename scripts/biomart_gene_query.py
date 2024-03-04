import argparse
from bioservices import BioMart

def get_gene_coordinates(gene_list_file, output_file, assembly='GRCh38'):
    server = BioMart()

    # Set the dataset and attributes
    dataset = 'hsapiens_gene_ensembl'
    server.add_dataset_to_xml(dataset)

    # Add attributes for canonical protein-coding exons and 5' UTR
    server.add_attribute_to_xml("ensembl_gene_id")
    server.add_attribute_to_xml("chromosome_name")
    server.add_attribute_to_xml("exon_chrom_start")
    server.add_attribute_to_xml("exon_chrom_end")
    server.add_attribute_to_xml("5_utr_start")
    server.add_attribute_to_xml("5_utr_end")

    # Read the list of genes from the input file
    with open(gene_list_file, 'r') as genes:
        for gene_symbol in genes:
            gene_symbol = gene_symbol.strip()

            # Add filter for the gene symbol
            server.add_filter_to_xml("hgnc_symbol", gene_symbol)

            # Execute the query
            xml_query = server.get_xml()
            result = server.query(xmlq=xml_query)

            # Process the result and append to the output file in BED format
            with open(output_file, 'a') as output:
                for line in result.split('\n'):
                    if line:
                        row = line.split('\t')
                        chromosome = row[1]
                        exon_start = row[2]
                        exon_end = row[3]
                        utr5_start = row[4]
                        utr5_end = row[5]

                        output.write(f"chr{chromosome}\t{exon_start}\t{exon_end}\n")

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Retrieve gene coordinates from BioMart.')
    parser.add_argument('--gene_list', required=True, help='Text file containing list of HGNC gene symbols, one per line')
    parser.add_argument('--output', required=True, help='Output BED file')
    args = parser.parse_args()

    # Run the function with provided arguments
    get_gene_coordinates(args.gene_list, args.output)
