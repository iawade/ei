import sys

bed_d = {}
gene_d = {}

script_dir = sys.path[0]

# May wish to redo this to accommodate your naming system, however first argument will also be fed
# into UKB.region_extraction.batch.sh as prefix for downstream analyses
with open(sys.argv[1] + '.UKB.GRCh38.regions.bed') as ipt:
    for line in ipt:
        line = line.strip().split('\t')
        chrom = int(line[0].replace('chr', '').replace('X', '23'))
        if chrom not in bed_d:
            bed_d[chrom] = []
        bed_d[chrom].append([int(line[1]), int(line[2]), line[3]])
        gene_d[line[-1]] = set()

for chrm in bed_d:
    bed_d[chrm] = sorted(bed_d[chrm], key=lambda x: (x[0], x[1]))

# UKB reference file listing pVCF block coordinates, (pvcf_blocks.txt)
# Again may need to adapt to your project/Snakemake structure
with open(script_dir + '/references/pvcf_blocks.txt') as ipt:
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
                    gene_d[exon[2]].add(file_name)

# Again may need to be adapted
with open(sys.argv[1] + '.UKB.files.txt', 'w') as opt:
    for gene in sorted(list(gene_d.keys())):
        if not gene_d[gene]:
            print('Coordinates could not be cross-referenced for', gene)
            opt.write(f'{gene}\t-\n')
        else:
            files = sorted(list(gene_d[gene]))
            for i in range(len(files)):
                opt.write(f'{gene}\t{i+1}\t{files[i]}\n')
