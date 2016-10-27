from data import parse_fna, parse_gbff
import re

END_CODONS = {'TAA', 'TAG', 'TGA'}
LOWER_CUTOFF = 50
UPPER_CUTOFF = 1400

def predict_by_length(seq, lower_bound=LOWER_CUTOFF, upper_bound=UPPER_CUTOFF):
    frames = [[] for i in range(3)]
    for offset in range(3):
        orf_start = 1 + offset
        for i in range(offset, len(seq), 3):
            if len(seq) - i >= 3: # check to make sure we don't run off the end of the genome
                codon = re.sub(r'[^AGTC]', 'T', seq[i:i+3].upper())
                if codon in END_CODONS:
                    frames[offset].append((orf_start, i+3))
                    orf_start = i+4
    return frames

if __name__ == '__main__':
    a = predict_by_length(parse_fna())
    over_1400 = 0
    under_50 = 0
    mid = 0
    for j in a:
        for i in j:
            if i[1] - i[0] > 1400:
                over_1400 += 1
            elif i[1] - i[0] < 50:
                under_50 += 1
            else:
                mid += 1
    print(sum([len(b) for b in a]))