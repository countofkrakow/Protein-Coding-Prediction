import re

def parse_fna(fname='GCF_000091665.1_ASM9166v1_genomic.fna'):
    gene = ''
    for line in open(fname, 'r'):
        if '>' in line:
            if gene:
                return gene
        else:
            gene += line.strip()
    return gene

def parse_gbff(fname='GCF_000091665.1_ASM9166v1_genomic.gbff'):
    f = open(fname, 'r')
    gene_pattern = re.compile("CDS\s+(\d+)\.\.(\d+)")
    gene_1 = f.read().split('ORIGIN')[0]
    return [(int(seq[0]), int(seq[1])) for seq in re.findall(gene_pattern, gene_1)]

# GCF_000091665.1_ASM9166v1_genomic.fna
# GCF_000091665.1_ASM9166v1_genomic.gbff

if __name__ == '__main__':
    p = parse_gbff('GCF_000091665.1_ASM9166v1_genomic.gbff')
    print(str([len(a) % 3 for a in p]))
