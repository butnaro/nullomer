import re,os,sys,glob
import mycustom 
from Bio.Seq import Seq 
from collections import Counter
from argparse import ArgumentParser
import pickle
import time

# Get sequences to run
def get_sequences(fasta):
	filed = fasta
	sequenceO = mycustom.FastaFile(filed)
	sequencesL = [ i.sequence.upper() for i in sequenceO ]
	sequencesL_rev_compl = []
	for i in sequencesL:
		seq = Seq(i)
		sequencesL_rev_compl+=[str(seq.reverse_complement())]
	return sequencesL, sequencesL_rev_compl

# Count motifs
def fast_counter(SequencesL,SequencesL_rev_compl, kmer_length, name): 
	Counts = Counter()
	for seq in SequencesL:
		for i in range(0, len(seq)+1-kmer_length):
			motif = seq[i:i+kmer_length]
			Counts[motif] += 1
	for seq2 in SequencesL_rev_compl:
		for im in range(0, len(seq2)+1-kmer_length):
			motif = seq2[im:im+kmer_length]
			Counts[motif] += 1
	most_common = Counts.most_common()
	with open(name+'.pickle', 'wb') as handle:
		pickle.dump(most_common, handle, protocol=pickle.HIGHEST_PROTOCOL)
	datafile=open(name + ".tsv","w")
	for i in most_common:
		datafile.write(str(i[0])+'\t'+str(i[1])+'\n')
	datafile.close()
	return

if __name__ == "__main__":
	parser = ArgumentParser()
	parser.add_argument("-f", "--fasta",
				   help="input fasta file")
	parser.add_argument("-k", "--kmer",
					nargs="*",
					type=int,
					default=[13,14,15],
					help="list of kmer lengths to count")
	args = parser.parse_args()
	(sequencesL, sequencesL_rev_compl) = get_sequences(args.fasta)
	for kmer_length in args.kmer:
		start = time.time()
		name = "test.fasta".split('.fasta')[0] + "_" + str(kmer_length) 
		fast_counter(sequencesL,sequencesL_rev_compl, kmer_length, name)
		end = time.time()
		print("Time elapsed for kmer of size " + str(kmer_length) + ": " + str(end - start))

