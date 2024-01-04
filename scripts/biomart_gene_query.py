import argparse
from bioservices import BioMart

def get_gene_coordinates(gene_symbol, output_file, assembly='GRCh38'):
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

    # Add filter for the gene symbol
    queries = gene_symbol
    server.add_filter_to_xml("hgnc_symbol", queries)

    # Execute the query
    xml_query = server.get_xml()
    result = server.query(xmlq=xml_query)

    # Print the entire result for investigation
    #print(result)

    # Process the result and write to the output file in BED format
    with open(output_file, 'w') as output:
        for line in result.split('\n'):
            if line:
                row = line.split('\t')
                gene_id = row[0]  
                chromosome = row[1]
                exon_start = row[2]
                exon_end = row[3]
                utr5_start = row[4]
                utr5_end = row[5]

                output.write(f"chr{chromosome}\t{exon_start}\t{exon_end}\n")
        #print(f"chr{chromosome}\t{utr5_start}\t{utr5_end}\t{gene_id}\t5_UTR")

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Retrieve gene coordinates from BioMart.')
    parser.add_argument('--gene', required=True, help='HGNC gene symbol')
    parser.add_argument('--output', required=True, help='Output BED file')
    args = parser.parse_args()

    # Run the function with provided arguments
    get_gene_coordinates(args.gene, args.output)
