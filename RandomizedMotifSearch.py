from random import randint
from GreedyMotifSearch import score,profile_most_probable_kmer

def profile_with_pseudocounts(motifs):
    prof = []
    for i in range(len(motifs[0])):
        col = ''.join([motifs[j][i] for j in xrange(len(motifs))])
        prof.append([float(col.count(nuc)+1)/float(len(col)+4) for nuc in 'ACGT'])
    return prof

def motifs_from_profile(profile, dna, k):
	return [profile_most_probable_kmer(seq,k,profile) for seq in dna]

def randomized_motif_search(dna,k,t):
	rand_ints = [randint(0,len(dna[0])-k) for a in range(t)]
	motifs = [dna_list[i][r:r+k] for i,r in enumerate(rand_ints)]

	best_score = [score(motifs), motifs]

	while True:
		current_profile = profile_with_pseudocounts(motifs)
		motifs = motifs_from_profile(current_profile, dna_list, k)
		current_score = score(motifs)
		if current_score < best_score[0]:
			best_score = [current_score, motifs]
		else:
			return best_score

if __name__ == '__main__':

	with open('randomized_motif') as input_data:
		k,t = map(int, input_data.readline().split())
		dna_list = [line.strip() for line in input_data.readlines()]

	best_motifs = [k*t, None]

	for repeat in range(1000):
		current_motifs = randomized_motif_search(dna_list,k,t)
		if current_motifs[0] < best_motifs[0]:
			best_motifs = current_motifs

	print len(best_motifs[1])
