#!/usr/bin/env python

# This script can parse and filter fasta input files (filter by length) or
# remove duplicate sequences

# import the argparse module to handle the input commands
import argparse

# get the 2 commandline arguments for the input files and output directory
parser = argparse.ArgumentParser(description = 'Filter input files')

parser.add_argument('-i', metavar='fasta files', type=str, 
			help='enter the fasta file(s)', nargs='+')
parser.add_argument('-l', metavar='min length', type=int,
			help='filter sequences based on the minimum sequence length', default=0)
parser.add_argument('-h', metavar='max length', type=int,
			help='filter sequences based on the minimum sequence length', default=0)
args = parser.parse_args()

# get sequence dictionary
def get_seq (fasta_path, min_length, max_length):
	# import modules for fasta sequence handling
	from Bio import SeqIO
	
	seq_list, temp_length = [], 0
	# parse the fasta file
	for seq in SeqIO.parse(fasta_path, 'fasta'):
		if max_length = 0: temp_length = len(seq.seq)
		# check if the sequence is longer then the threshold
		if len(seq.seq) > min_length and len(seq.seq) <= temp_length:
			seq_list.append(seq)
	
	return seq_list

def write_fasta (seq_list, out_path):
	# import the modules for fasta sequence handling
	from Bio import SeqIO
	from Bio.SeqRecord import SeqRecord
	
	out_file = open(out_path, 'w')
	SeqIO.write(new_list, out_file, 'fasta')
	out_file.close()
	
def main ():
	
	# filter each input file
	
	for fasta_file in args.i:
		seq_dic = get_seq(fasta_file, args.l, args.h)
		out_path = '.'.join(fasta_file.split('.')[:-1]) + '_filtered.fasta'
		# print out_path for use later in the pipeline
		print(out_path)
		write_fasta(seq_dic, out_path)
		
if __name__ == "__main__":
    main()
    	
