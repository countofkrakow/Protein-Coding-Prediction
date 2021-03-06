from data import parse_fna, parse_gbff
from math import log2
import re

END_CODONS = {'TAA', 'TAG', 'TGA'}
LOWER_CUTOFF = 50
UPPER_CUTOFF = 1400

def predict_by_length(seq, lower=LOWER_CUTOFF, upper=UPPER_CUTOFF):
    predictions = [[[] for j in range(3)] for i in range(3)]
    for offset in range(3):
        orf_start = 1 + offset
        for i in range(offset, len(seq), 3):
            if len(seq) - i >= 3: # check to make sure we don't run off the end of the genome
                codon = re.sub(r'[^AGTC]', 'T', seq[i:i+3].upper())
                if codon in END_CODONS:
                    start_pos = orf_start
                    end_pos = i + 3
                    if end_pos - start_pos < lower:
                        predictions[0][offset].append((start_pos, end_pos))
                    elif end_pos - start_pos > upper:
                        predictions[2][offset].append((start_pos, end_pos))
                    else:
                        predictions[1][offset].append((start_pos, end_pos))
                    orf_start = i+4
    return predictions

def predict_by_threshold(seq, threshold=50):
    predictions = [[] for j in range(2)]
    for offset in range(3):
        orf_start = 1 + offset
        for i in range(offset, len(seq), 3):
            if len(seq) - i >= 3: # check to make sure we don't run off the end of the genome
                codon = re.sub(r'[^AGTC]', 'T', seq[i:i+3].upper())
                if codon in END_CODONS:
                    start_pos = orf_start
                    end_pos = i + 3
                    if end_pos - start_pos < threshold:
                        predictions[0].append((start_pos, end_pos))
                    else:
                        predictions[1].append((start_pos, end_pos))
                    orf_start = i+4
    return predictions

def check_predictions(predictions, actual):
    # p = set([orf[1] for offsets in predictions for ranges in offsets for orf in ranges])
    p = set([range[1] for range in predictions])
    true_positives = len(p & actual)
    false_positives = len(p - actual)
    false_negatives = len(actual - p)
    return true_positives, false_positives, false_negatives

def train_markov_counts(sequences, chain_len=3):
    counts = {}
    for seq in sequences:
        seq = seq[:-3]
        for i in range(chain_len, len(seq)):
            markov = re.sub(r'[^AGTC]', 'T', seq[i - chain_len:i].upper())
            unk = re.sub(r'[^AGTC]', 'T', seq[i].upper())
            if markov not in counts:
                counts[markov] = {}
            if unk not in counts[markov]:
                counts[markov][unk] = 0
            counts[markov][unk] += 1
    return counts

def train_markov_model(valid_orfs, invalid_orfs):
    P = train_markov_counts(valid_orfs)
    Q = train_markov_counts(invalid_orfs)

    map = {0:'A', 1:'C', 2:'G', 3:'T'}
    tableP = [[0]*4 for i in range(4)]
    tableQ = [[0]*4 for i in range(4)]

    for x in range(4):
        for y in range(4):
            markov = 'A' + map[x] + map[y]
            tableP[x][y] = (P[markov]['A'], sum([P[markov][a] for a in P[markov]]))
            tableQ[x][y] = (Q[markov]['A'], sum([Q[markov][a] for a in Q[markov]]))
    return P, Q, tableP, tableQ

def print_markov_table(counts):
    map = {0:'A', 1:'C', 2:'G', 3:'T'}
    print('x/y  {0:14}{1:14}{2:14}{3:14}'.format('A', 'C', 'G', 'T'))
    for i, row in enumerate(counts):
        print(map[i] + ' '*4 + '{0:14}{1:14}{2:14}{3:14}'.format(str(row[0]), str(row[1]), str(row[2]), str(row[3])))

def markov_score(P, Q, seq, chain_len=3):
    mm = 0
    for j in range(chain_len, len(seq)):
        markov = seq[j-chain_len:j]
        nuc = seq[j]
        pp = log2(P[markov][nuc]) - log2(sum([P[markov][i] for i in P[markov]]))
        qq = log2(Q[markov][nuc]) - log2(sum([Q[markov][i] for i in Q[markov]]))
        mm += pp - qq
    return mm

if __name__ == '__main__':
    raw_data = parse_fna()
    genes = parse_gbff()
    predictions = predict_by_length(raw_data)

    # Make predictions based on length
    print('1.a ORF\'s per reading frame:\n\t0: %d\n\t1: %d\n\t2: %d' %  (
        sum([len(offsets[0]) for offsets in predictions]),
        sum([len(offsets[1]) for offsets in predictions]),
        sum([len(offsets[2]) for offsets in predictions])))
    print('1.b ORF\'s shorter than 50 nucleotides: ' + str(sum([len(orf) for orf in predictions[0]])))
    print('1.c ORF\'s longer than 1400 nucleotides: ' + str(sum([len(orf) for orf in predictions[2]])))

    # check predictions made based on length
    all_pred = [r for length in predictions for ranges in length for r in ranges]

    actual = set([g[1] for g in genes])
    tru_pos, false_pos, false_neg = check_predictions(all_pred, actual)
    print('1.d CDS\'s found in GenBank: ' + str(tru_pos))

    # Markov count tables
    valid_orfs = [raw_data[start-1:end] for offset in predictions[2] for start, end in offset]
    invalid_orfs = [raw_data[start-1:end] for offset in predictions[0] for start, end in offset]
    P, Q, tableP, tableQ = train_markov_model(valid_orfs, invalid_orfs)
    print('1.e Markov chain counts')
    print('P counts:')
    print_markov_table(tableP)
    print('Q counts:')
    print_markov_table(tableQ)

    lowest_short_orfs = []
    lowest_long_orfs = []
    for offset in predictions[0]:
        lowest_short_orfs += offset
    lowest_short_orfs.sort()
    for offset in predictions[2]:
        lowest_long_orfs += offset
    lowest_long_orfs.sort()

    # ORF's lower than 50
    print('SUMMARY FOR FIRST 5 ORFS OF LENGTH < 50')
    print('#     Start      Length     MM Score           CDS?')
    for i in range(5):
        start = lowest_short_orfs[i][0] - 1
        end = lowest_short_orfs[i][1]
        length = end - start
        seq = raw_data[start:end]
        mm = markov_score(P, Q, seq)
        print('{0:d}{1:6d}{2:12d}{3:22.11f}      {4:18}'.format(i+1, start, length, mm, str(end in actual)))

    # ORF's higher than 1400
    print('SUMMARY FOR FIRST 5 ORFS OF LENGTH > 1400')
    print('#     Start      Length     MM Score           CDS?')
    for i in range(5):
        start = lowest_long_orfs[i][0] - 1
        end = lowest_long_orfs[i][1]
        length = end - start
        seq = raw_data[start:end]
        mm = markov_score(P, Q, seq)
        print('{0:d}    {1:6d}{2:10d}{3:21.11f}     {4:18}'.format(i+1, start, length, mm, str(end in actual)))

