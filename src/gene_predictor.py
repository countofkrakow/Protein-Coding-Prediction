from data import parse_fna, parse_gbff
import re

END_CODONS = {'TAA', 'TAG', 'TGA'}
LOWER_CUTOFF = 50
UPPER_CUTOFF = 1400

def predict_by_length(seq, lower=LOWER_CUTOFF, upper=UPPER_CUTOFF):
    frames = [[[] for j in range(3)] for i in range(3)]
    for offset in range(3):
        orf_start = 1 + offset
        for i in range(offset, len(seq), 3):
            if len(seq) - i >= 3: # check to make sure we don't run off the end of the genome
                codon = re.sub(r'[^AGTC]', 'T', seq[i:i+3].upper())
                if codon in END_CODONS:
                    start_pos = orf_start
                    end_pos = i + 3
                    if end_pos - start_pos < lower:
                        frames[0][offset].append((start_pos, end_pos))
                    elif end_pos - start_pos > upper:
                        frames[2][offset].append((start_pos, end_pos))
                    else:
                        frames[1][offset].append((start_pos, end_pos))
                    orf_start = i+4
    return frames

def check_predictions(predictions, actual):
    p = set([orf[1] for offsets in predictions for ranges in offsets for orf in ranges])
    a = set([a[1] for a in actual])
    true_positives = len(p & a)
    false_positives = len(p - a)
    false_negatives = len(a - p)
    return true_positives, false_positives, false_negatives

def train_markov_model(data, actual, chain_len=3):
    # mapping: 'ATC' -> {
    P = {}
    Q = {}


if __name__ == '__main__':
    fna = parse_fna()
    actual = parse_gbff()
    predictions = predict_by_length(fna)
    a = [orf[1] for range in predictions[0] for orf in range]

    print('1.a ORF\'s per reading frame:\n\t0: %d\n\t1: %d\n\t2: %d' %  (\
        sum([len(offsets[0]) for offsets in predictions]),
        sum([len(offsets[1]) for offsets in predictions]),
        sum([len(offsets[2]) for offsets in predictions])))


    print('1.b ORF\'s shorter than 50 nucleotides: ' + str(sum([len(orf) for orf in predictions[0]])))
    print('1.c ORF\'s longer than 1400 nucleotides: ' + str(sum([len(orf) for orf in predictions[2]])))

    tru_pos, false_pos, false_neg = check_predictions(predictions, actual)
    print('1.d CDS\'s found in GenBank: ' + str(tru_pos))
    print(check_predictions(predictions, actual))