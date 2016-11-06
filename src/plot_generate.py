import gene_predictor
from data import parse_fna, parse_gbff
import re
from statistics import median

# generates ROC curve based off orf length
def generate_roc_length(actual, raw_data):
    fp = []
    tp = []
    end_pos = []
    for i in range(50, 1410, 10):
        preds = gene_predictor.predict_by_threshold(raw_data, i)
        true_positives, false_positives, false_negatives = gene_predictor.check_predictions(preds[1], actual)

        tru_pos_rate = true_positives / len(actual)
        false_pos_rate = false_positives / (len(preds[0]) + len(preds[1]) - len(actual))

        tp.append(tru_pos_rate)
        fp.append(false_pos_rate)
    return tp, fp

# generates ROC curve based off markov scores
def generate_roc_markov(actual, raw_data, P, Q):
    fp = []
    tp = []
    all_orfs = gene_predictor.predict_by_threshold(raw_data)
    all_orfs = all_orfs[0] + all_orfs[1]
    for markov_score in range(100):
        preds = []
        total = 0
        for r in all_orfs:
            start = r[0] - 1
            end = r[1]
            seq = raw_data[start:end]
            mm = gene_predictor.markov_score(P, Q, seq)
            if mm > (float(markov_score)/5):
                preds.append(r)
            total += 1
        true_positives, false_positives, false_negatives = gene_predictor.check_predictions(preds, actual)

        tru_pos_rate = true_positives / len(actual)
        false_pos_rate = false_positives / (total - len(actual))

        tp.append(tru_pos_rate)
        fp.append(false_pos_rate)
    return tp, fp

# makes predictions off the dr flashbulb line
def calc_line(data):
    # Read files
    with open('valid.csv', 'r') as f:
        valid = []
        for line in f:
            line = line.split(',')
            valid.append((float(line[0]), int(line[1])))

    with open('invalid.csv', 'r') as f:
        invalid = []
        for line in f:
            line = line.split(',')
            invalid.append((float(line[0]), int(line[1])))

    with open('unk.csv', 'r') as f:
        unk = []
        for line in f:
            line = line.split(',')
            unk.append((float(line[0]), int(line[1])))

    # calculate the two data
    Ay = median([p[0] for p in invalid])
    Ax = median([p[1] for p in invalid])
    By = median([p[0] for p in valid])
    Bx = median([p[1] for p in valid])
            

    # line = mx + b
    m = (By - Ay) / (Bx - Ax)

    # b = Ay - m*Ax
    threshold_m = -1/m
    with open('length_vs_markov_threshold.txt', 'w') as f:

        for i in range(100):
            threshold = i/100
            threshold_x = Ax + threshold*(Bx - Ax)
            threshold_y = Ay + threshold*(By - Ay)
            threshold_b = threshold_y - threshold_m*threshold_x
            true_positives, false_positives, true_negatives, false_negatives = predict_flashbulb(threshold_b, threshold_m, data)
            fpr = false_positives / (true_negatives + false_positives)
            tpr = true_positives / (true_positives + false_negatives)
            f.write(str(false_positives) + ',' + str(fpr) + ',' + str(tpr) + ',' + str(i) + '\n')


# Predicts for one given flashbulb line
def predict_flashbulb(b, m, data):
    true_positives = 0
    false_positives = 0
    false_negatives = 0
    true_negatives = 0

    for x, y, pos in data:
        # pos prediction
        if y - m*x >= b:
            if pos:
                true_positives += 1
            else:
                false_positives += 1
        # neg prediction
        else:
            if pos:
                false_negatives += 1
            else:
                true_negatives += 1

    return true_positives, false_positives, true_negatives, false_negatives

# writes down length and markov score pairs for each orf
def graphVals(valid_orfs, invalid_orfs, unknown_orfs, actual, raw_data, P, Q):
    # for plotting
    markov_scores = [[], [], []]
    orf_scores = [[], [], []]
    positions = [[], [], []]
    for orf in invalid_orfs:
        score = gene_predictor.markov_score(P, Q, orf)
        markov_scores[0].append(score)
        orf_scores[0].append(len(orf))
        positions[0].append(orf)

    for orf in valid_orfs:
        score = gene_predictor.markov_score(P, Q, orf)
        markov_scores[2].append(score)
        orf_scores[2].append(len(orf))
        positions[2].append(orf)

    for orf in unknown_orfs:
        score = gene_predictor.markov_score(P, Q, orf)
        markov_scores[1].append(score)
        orf_scores[1].append(len(orf))
        positions[1].append(orf)

    with open('unk.csv', 'w') as f:
        for i, score in enumerate(markov_scores[1]):
            f.write(str(score) + ',' + str(orf_scores[1][i]) + '\n')

    with open('valid.csv', 'w') as f:
        for i, score in enumerate(markov_scores[2]):
            f.write(str(score) + ',' + str(orf_scores[2][i]) + '\n')

    with open('invalid.csv', 'w') as f:
        for i, score in enumerate(markov_scores[0]):
            f.write(str(score) + ',' + str(orf_scores[0][i]) + '\n')

# calculates markov score, orf length, and sequence end for each data point
def calc_markov_length_end(raw_data, P, Q, actual):
    pos_data = gene_predictor.predict_by_threshold(raw_data)
    pos_data = pos_data[0] + pos_data[1]
    v = []
    for orf in pos_data:
        start = orf[0]
        end = orf[1]
        seq = raw_data[start-1:end]
        length = len(seq)
        markov = gene_predictor.markov_score(P, Q, seq)
        is_valid = end in actual
        v.append((length, markov, is_valid))
    return v


if __name__ == '__main__':
    raw_data = re.sub(r'[^AGTC]', 'T', parse_fna().upper())
    genes = parse_gbff()
    predictions = gene_predictor.predict_by_length(raw_data)
    valid_orfs = [raw_data[start-1:end] for offset in predictions[2] for start, end in offset]
    invalid_orfs = [raw_data[start-1:end] for offset in predictions[0] for start, end in offset]
    unk_orfs = [raw_data[start-1:end] for offset in predictions[1] for start, end in offset]
    P, Q, tableP, tableQ = gene_predictor.train_markov_model(valid_orfs, invalid_orfs)

    # flags
    write_markov_roc_csv = False
    write_markov_vs_orf = False
    write_length_csv = False
    write_markov_csv = False
    uh = True

    actual = set([a[1] for a in genes])

    if write_markov_vs_orf:
        tp_l, fp_l = generate_roc_length(actual, raw_data)
        tp_m, fp_m = generate_roc_markov(actual, raw_data, P, Q)
        with open("markov_roc_curve.csv", 'w') as f:
            for i in range(len(tp_m)):
                f.write(str(tp_m[i]) + "," + str(fp_m[i]) + '\n')
        with open("length_roc_curve.csv", 'w') as f:
            for i in range(len(tp_l)):
                f.write(str(tp_l[i]) + "," + str(fp_l[i]) + '\n')

    if write_markov_roc_csv:
        graphVals(valid_orfs, invalid_orfs, unk_orfs, actual, raw_data, P, Q)

    if uh:
        v = calc_markov_length_end(raw_data, P, Q, actual)
        calc_line(v)

if write_length_csv:
    res = []
    for i in range(50, 1500, 10):
        pred = gene_predictor.predict_by_threshold(raw_data, i)
        tru_pos, false_pos, false_neg = gene_predictor.check_predictions(pred[1], actual)
        res.append((str(i), str(float(tru_pos)/(tru_pos+false_neg)), str(tru_pos), str(false_pos)))
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
        res.append((str(float(markov_score)/5), str(float(tru_pos)/(tru_pos + false_neg)), str(tru_pos), str(false_pos)))
    f = open('markov.csv', 'w')
    for point in res:
        f.write(point[0] + ', ' + point[1] + ', ' + point[2] + ', ' + point[3] + '\n')
    f.close()

