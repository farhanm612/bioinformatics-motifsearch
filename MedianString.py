from operator import ne
from itertools import product


def hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError('Undefined for sequences of unequal length.')
    return sum(map(ne, seq1, seq2))


def motif_score(pattern, motif):
	return min([hamming_distance(motif[i:i+len(pattern)], pattern) for i in range(len(motif)-len(pattern)+1)])

if __name__ == "__main__":
    with open('median_string.txt') as input_data:
        k = int(input_data.readline())
        dna_list = [line.strip() for line in input_data.readlines()]

    best_pattern = [k * len(dna_list) + 1]
    for pattern in product('ACGT', repeat=k):
        current_score = sum([motif_score(''.join(pattern), dna) for dna in dna_list])
        if current_score < best_pattern[0]:
            best_pattern = [current_score, ''.join(pattern)]

    print best_pattern[1]
