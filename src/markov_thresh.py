import gene_predictor
from data import parse_fna, parse_gbff
import re
import matplotlib.pyplot as plt

raw_data = re.sub(r'[^AGTC]', 'T', parse_fna().upper())
genes = parse_gbff()
predictions = gene_predictor.predict_by_length(raw_data)
valid_orfs = [raw_data[start-1:end] for offset in predictions[2] for start, end in offset]
invalid_orfs = [raw_data[start-1:end] for offset in predictions[0] for start, end in offset]
P, Q, tableP, tableQ = gene_predictor.train_markov_model(valid_orfs, invalid_orfs)

# flags

write_markov_roc_csv = False
write_markov_vs_orf = True

actual = set([a[1] for a in genes])

res = []
for markov_score in range(100):
    all_orfs = gene_predictor.predict_by_threshold(raw_data)
    all_orfs = all_orfs[0] + all_orfs[1]
    pred = []
    for range in all_orfs:
        start = range[0] - 1
        end = range[1]
        seq = raw_data[start:end]
        mm = gene_predictor.markov_score(P, Q, seq)
        if mm > (float(markov_score)/5):
            pred.append(range)

    tru_pos, false_pos, false_neg = gene_predictor.check_predictions(pred, actual)
    res.append((str(float(markov_score)/5), str(tru_pos), str(false_pos), str(float(tru_pos)/(tru_pos + false_pos))))
f = open('markov.csv', 'w')
for point in res:
    f.write(point[0] + ', ' + point[1] + ', ' + point[2] + ', ' + point[3] + '\n')
f.close()