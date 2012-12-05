# The pick_otu_rep.py script creates a representative sequence for a cluster of
# sequences. It can pick a sequence based on three different methods: 'random',
# a random sequence will be picked as the representative sequence, 'consensus',
# the representative sequence is based on the consensus sequence for the cluster
# and 'combined' wich takes an x number of random sequences from a cluster and
# creates a consensus sequence for these ('combined' helps to deal with large
# clusters, since a full sequence alignment will greatly slow the script down).


# import the argparse module to handle the input commands
import argparse
import sys

# get the commandline arguments for the input files and output directory and user settings
parser = argparse.ArgumentParser(description = 'retrieve a represented otu from each cluster')

parser.add_argument('-i', metavar='fasta file(s)', type=str, 
			help='enter the fasta files that where clustered', nargs='+')
parser.add_argument('-o', metavar='output file', type=str, 
			help='enter the output file')
parser.add_argument('-c', metavar='cluster file', type=str, 
			help='enter the cluster (otu) file')
parser.add_argument('-m', metavar='minimal cluster size', type=int,
			help='enter the minimal cluster size')
parser.add_argument('-s', metavar='select method', type=str, 
			help='the method used to select a representative sequence for a cluster: \'random\' (default), \'consensus\' or \'combined\', please note: using the consensus option on large cluster can take a lot of time and computer power to finnish, use \'combined\' when more speed is needed', default='random')
parser.add_argument('-r', metavar='number of random sequences used for the consensus', type=int,
			help='the number of random sequences that will be used to create a consensus sequence (select method: \'combined\' only', default=10)
parser.add_argument('-p', metavar='cluster program  used', type=str,
			help='the program used for the cluster analysis')
args = parser.parse_args()

def extract_otu (otufile):
	# retrieve the otu clusters that meet the size restriction
	
	otu_seq_dic, sequence_header_dic = {}, {}

	# parse through the OTU file
	for line in open(otufile, 'r'):
		# split each OTU into a set of fasta header that make up
		# the OTU, if the number of header matches the size threshold
		# the cluster is stored
		array = line.replace('\n', '').split('\t')
		#if len(array[1:]) >= minsize: 
		otu_seq_dic[array[0]] = array[1:]
		for header in array[1:]:
			sequence_header_dic[header] = [array[0], len(array[1:])]

	return [otu_seq_dic, sequence_header_dic]

def extract_seq (seq_file):
	# The SeqIO module is imported from Biopython to ease the parsing of fasta files
	from Bio import SeqIO
	
	seq_dic = {}

	# parse through a list of fasta file
	for fasta_file in seq_file:
		# Create a dictionary of fasta sequences with the headers as keys
		for seq_record in SeqIO.parse(fasta_file, 'fasta'):
			seq_dic[seq_record.id] = seq_record.seq
	
	# return the dictonary containing the sequences from the list of input files
	return seq_dic

#def extract_consensus (seq_file, otu_seq_dic, program):
#	# the SeqIO module is imported from Biopython to parse the fasta file
#	from Bio import SeqIO
#
#	seq_dic, header = {}, ''
#
#	# parse through the consensus file
#	for seq_record in SeqIO.parse(seq_file[0], 'fasta'):
#		if program == 'usearch': header = seq_record.id.replace('Cluster')
#		elif program == 'usearch_old': header = seq_record.id.replace('Cluster')

def write_results (sequence, header, cluster, clust_length, out_path):
	# import a set of Biopython modules for sequence handling and writing
	from Bio import SeqIO
	from Bio.SeqRecord import SeqRecord	

	# write a fasta sequence to the output file, the header of the sequence
	# contains information about the cluster it is based on and the length of the cluster
	out_file = open(out_path, 'a')
	seq = SeqRecord(sequence, id=(header + '_cluster_#:' + str(cluster) + '_length_cluster:' + str(clust_length) + '_'), description='')
	SeqIO.write(seq, out_file, 'fasta')
	out_file.close()

def get_rand_seq (seq_dic, otu_seq_dic, out_path, min_size):
	from random import choice

	# keys = otu_seq_dic.keys()
	
	# parse through the dictionary containing the headers for each otu
	for item in otu_seq_dic[0]: #keys:
		# randomly select 1 header to represent the otu
		header = choice(otu_seq_dic[0][item])
		# the header + sequence are writen to the output fasta file
		if len(otu_seq_dic[0][item]) >= min_size:
			write_results(seq_dic[header], header, str(item), len(otu_seq_dic[0][item]), out_path)

def get_cons_seq (seq_dic, otu_seq_dic, out_path, program, min_size):
	# add aditional information to the cluster consensus files (cluster number and size)

	header, length = '', 0

	# parse through the consensus sequence dictionary	
	for seq in seq_dic:
		if program == 'octupus': header = seq.replace('OCTU','')
		elif program == 'usearch_old': header = seq.replace('Cluster','')
		elif program == 'usearch': header = otu_seq_dic[1][seq.replace('centroid=','').split(';')[0]][0]
		# write the sequence with some aditional information on the cluster size
		if len(otu_seq_dic[0][header]) >= min_size:
			write_results(seq_dic[seq], seq, header, len(otu_seq_dic[0][header]), out_path)
	

	#for item in otu_seq_dic:
	#	header = item
	#	# write the sequence with some aditional information on the cluster size
	#	write_results(seq_dic[header], header, str(item), len(otu_seq_dic[item]), out_path)

	
def main ():
	
	# check which method needs to be used to retrieve the representative sequence
	if args.s == 'random':
		get_rand_seq(extract_seq(args.i), extract_otu(args.c), args.o, args.m)
	elif args.s == 'consensus':
		get_cons_seq(extract_seq(args.i), extract_otu(args.c), args.o, args.p, args.m)	

if __name__ == "__main__":
    main()


