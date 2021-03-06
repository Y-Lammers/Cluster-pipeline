#!/usr/bin/env python

# Pick a representative sequence for the clusters, this
# can either be a consensus sequence (if supported by the
# cluster program) or a random sequence from the cluster


# import the argparse module to handle the input commands
import argparse
import sys

# get the commandline arguments for the input files and output directory and user settings
parser = argparse.ArgumentParser(description = 'retrieve a represented otu from each cluster')

parser.add_argument('-i', metavar='fasta file', type=str, 
			help='enter the clustered fasta file')
parser.add_argument('-o', metavar='output file', type=str, 
			help='enter the output file')
parser.add_argument('-c', metavar='cluster file', type=str, 
			help='enter the cluster (otu) file')
parser.add_argument('-s', '--min_size', metavar='minimal cluster size', dest='s', type=int,
			help='enter the minimal cluster size')
parser.add_argument('-m','--method', metavar='consensus / random', dest='m', type=str, 
			help='the method used to select a representative sequence for a cluster: consensus or random (default consensus) * note the consensus option is not available for cd-hit or tgicl, it will automatically change to the random setting', default='consensus')
parser.add_argument('-p', metavar='cluster program  used', type=str,
			help='the program used for the cluster analysis (default: octupus)', default='octupus')
args = parser.parse_args()

def extract_otu ():
	# retrieve the otu clusters that meet the size restriction
	
	otu_seq_dic, sequence_header_dic = {}, {}

	# parse through the OTU file
	for line in open(args.c, 'r'):
		# split each OTU into a set of fasta header that make up
		# the OTU, if the number of header matches the size threshold
		# the cluster is stored
		array = line.strip().split('\t')
		otu_seq_dic[array[0]] = array[1:]
		for header in array[1:]:
			sequence_header_dic[header] = [array[0], len(array[1:])]

	return [otu_seq_dic, sequence_header_dic]

def extract_seq ():
	# The SeqIO module is imported from Biopython to ease the parsing of fasta files
	from Bio import SeqIO
	
	seq_dic = {}

	# Create a dictionary of fasta sequences with the headers as keys
	for seq_record in SeqIO.parse(args.i, 'fasta'):
		seq_dic[seq_record.id] = seq_record.seq
	
	# return the dictonary containing the sequences from the list of input files
	return seq_dic

def write_results (sequence, header, cluster, clust_length):
	# import a set of Biopython modules for sequence handling and writing
	from Bio import SeqIO
	from Bio.SeqRecord import SeqRecord	

	# write a fasta sequence to the output file, the header of the sequence
	# contains information about the cluster it is based on and the length of the cluster
	out_file = open(args.o, 'a')
	seq = SeqRecord(sequence, id=(header + '_cluster_#:' + str(cluster) + '_length_cluster:' + str(clust_length) + '_'), description='')
	SeqIO.write(seq, out_file, 'fasta')
	out_file.close()


def get_rand_seq (seq_dic, otu_seq_dic):
	from random import choice

	# keys = otu_seq_dic.keys()
	
	# parse through the dictionary containing the headers for each otu
	for item in otu_seq_dic[0]: #keys:
		# randomly select 1 header to represent the otu
		header = choice(otu_seq_dic[0][item])
		# the header + sequence are writen to the output fasta file
		if len(otu_seq_dic[0][item]) >= args.s:
			write_results(seq_dic[header], header, str(item), len(otu_seq_dic[0][item]))


def get_cons_seq (seq_dic, otu_seq_dic):
	# add aditional information to the cluster consensus files (cluster number and size)

	header, length = '', 0

	# parse through the consensus sequence dictionary	
	for seq in seq_dic:

		# based on the used cluster program, retrieve the correct cluster number from the 
		# consensus file
		if args.p == 'octupus': header = seq.replace('OCTU','')
		elif args.p == 'usearch_old': 
			header = str(int(seq.replace('Cluster',''))+1)
		elif args.p == 'usearch': 
			header = str(int(otu_seq_dic[1][seq.replace('centroid=','').split(';')[0]][0])+1)
		
		# write the sequence with some aditional information on the cluster size
		if len(otu_seq_dic[0][header]) >= args.s:
			write_results(seq_dic[seq], seq, header, len(otu_seq_dic[0][header]))


def main ():

	# check which method needs to be used to retrieve the representative sequence
	if args.m == 'random':
		get_rand_seq(extract_seq(), extract_otu())
	elif args.m == 'consensus':
		get_cons_seq(extract_seq(), extract_otu())

if __name__ == "__main__":
	main()


