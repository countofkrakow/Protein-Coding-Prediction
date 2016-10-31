import gene_predictor
from data import parse_fna, parse_gbff
import re

raw_data = re.sub(r'[^AGTC]', 'T', parse_fna().upper())
genes = parse_gbff()
predictions = gene_predictor.predict_by_length(raw_data)
valid_orfs = [raw_data[start-1:end] for offset in predictions[2] for start, end in offset]
invalid_orfs = [raw_data[start-1:end] for offset in predictions[0] for start, end in offset]
P, Q, tableP, tableQ = gene_predictor.train_markov_model(valid_orfs, invalid_orfs)

write_length_csv = False
write_markov_csv = False
write_markov_roc_csv = True
write_length_roc_csv = True
actual = set([a[1] for a in genes])

if write_markov_roc_csv:
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
        res.append((str(float(false_pos)/(tru_pos + false_pos)), str(float(tru_pos)/(tru_pos + false_pos))))
    f = open('markov_roc.csv', 'w')
    for point in res:
        f.write(point[0] + ',' + point[1] + '\n')
    f.close()

if write_length_roc_csv:
    res = []
    for i in range(50, 1500, 10):
        pred = gene_predictor.predict_by_threshold(raw_data, i)
        tru_pos, false_pos, false_neg = gene_predictor.check_predictions(pred[1], actual)
        res.append((str(float(false_pos)/(tru_pos + false_pos)), str(float(tru_pos)/(tru_pos + false_pos))))
        #print(str(tru_pos) + " " + str(false_pos) + " " + str(false_neg))

    f = open('length_roc.csv', 'w')
    for point in res:
        f.write(point[0] + ',' + point[1] + '\n')
    f.close()

# Create csv file of length prediction
if write_length_csv:
    res = []
    for i in range(50, 1500, 10):
        pred = gene_predictor.predict_by_threshold(raw_data, i)
        tru_pos, false_pos, false_neg = gene_predictor.check_predictions(pred[1], actual)
        res.append((str(i), str(float(tru_pos)/(tru_pos+false_pos)), str(tru_pos), str(false_pos)))
        #print(str(tru_pos) + " " + str(false_pos) + " " + str(false_neg))

    f = open('length.csv', 'w')
    for point in res:
        f.write(point[0] + ', ' + point[1] + ', ' + point[2] + ', ' + point[3] + '\n')
    f.close()

if write_markov_csv:
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




# Create csv file of markov prediction
#f = open('markov.csv', 'w')
