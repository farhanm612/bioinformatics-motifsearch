from random import randint
from GreedyMotifSearch import score,profile_most_probable_kmer
from RandomizedMotifSearch import profile_with_pseudocounts

def gibbs_sampler(dna,k,t,N):
	rand_ints = [randint(0,len(dna[0])-k) for a in xrange(t)]
	motifs = [dna_list[i][r:r+k] for i,r in enumerate(rand_ints)]
	best_score = [score(motifs), motifs]
	for i in range(N):
		r = randint(0,t-1)
		current_profile = profile_with_pseudocounts([motif for index, motif in enumerate(motifs) if index!=r])
		motifs = [profile_most_probable_kmer(dna[index],k,current_profile) if index == r else motif for index,motif in enumerate(motifs)]
		current_score = score(motifs)
		if current_score < best_score[0]:
			best_score = [current_score, motifs]

	return best_score


if __name__ == '__main__':

	with open('Gibbs') as input_data:
		k,t,N = map(int, input_data.readline().split())
		dna_list = [line.strip() for line in input_data.readlines()]
	best_motifs = [k*t, None]

	for repeat in range(20):
		current_motifs = gibbs_sampler(dna_list,k,t,N)
		if current_motifs[0] < best_motifs[0]:
			best_motifs = current_motifs

	print best_motifs[1]