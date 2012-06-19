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
args = parser.parse_args()

def size_filt (otufile, minsize):
	# retrieve the otu clusters that meet the size restriction
	
	otu_seq_dic = {}	

	# parse through the OTU file
	for line in open(otufile, 'r'):
		# split each OTU into a set of fasta header that make up
		# the OTU, if the number of header matches the size threshold
		# the cluster is stored
		array = line.replace('\n', '').split('\t')
		if len(array[1:]) >= minsize: 
			otu_seq_dic[array[0]] = array[1:]

	return otu_seq_dic

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

def find_muscle_path ():
	# find the path to the muscle program
	path_file = open('/'.join(sys.argv[0].split('/')[:-1])+'/paths.txt', 'r')
	muscle_path = [line.split('\t')[1] for line in path_file if 'muscle' in line][0].replace('\n','')
		
	return muscle_path

def clean_up (file_path):
	from subprocess import call
	
	# removes a temporary file
	p = call(['rm', file_path])

def get_consensus ():
	# Imports the necessairy biopython modules to handle multiple sequence
	# alignment files and fasta files
	from Bio import AlignIO
	from Bio.Align import AlignInfo	
	#from Bio.SeqRecord import SeqRecord	

	# read the multiple sequence alignment file
	alignment = AlignInfo.SummaryInfo(AlignIO.read('temp_align.txt', 'fasta'))
	# return a consensus based on the alignment
	sequence = alignment.dumb_consensus(ambiguous='N')
	return sequence
	
def get_rand_seq (seq_dic, otu_seq_dic, out_path):
	from random import choice

	# keys = otu_seq_dic.keys()
	
	# parse through the dictionary containing the headers for each otu
	for item in otu_seq_dic: #keys:
		# randomly select 1 header to represent the otu
		header = choice(otu_seq_dic[item])
		# the header + sequence are writen to the output fasta file
		write_results(seq_dic[header], header, str(item), len(otu_seq_dic[item]), out_path)

def get_cons_seq (seq_dic, otu_seq_dic, seq_number, out_path):
	from subprocess import call
	from random import sample
	
	# get the headers that are present in each OTU and the path to the muscle executable
	keys, muscle_path = otu_seq_dic.keys(), find_muscle_path()
	
	# parse through the dictionary containing the headers for each otu
	for item in keys:	
		# write either all sequences or a subset to the 'temp_in.fasta' file
		if len(otu_seq_dic[item]) > seq_number and seq_number != False:		
			for header in sample(otu_seq_dic[item], seq_number):
				write_results(seq_dic[header], header, 1, len(otu_seq_dic[item]), 'temp_in.fasta')
		else:
			for header in otu_seq_dic[item]:
				write_results(seq_dic[header], header, 1, len(otu_seq_dic[item]), 'temp_in.fasta')
				
		# align the newly created 'temp_in.fasta' file with muscle
		p = call([muscle_path, '-in', 'temp_in.fasta', '-out', 'temp_align.txt', '-quiet'])
		# remove the 'temp_in.fasta' file since this file is no longer needed
		clean_up('temp_in.fasta')

		# Get the consensus sequence for the multiple sequence alignment, and write it to
		# the output file.
		write_results(get_consensus(), header, str(item), len(otu_seq_dic[item]), out_path)
		# remove the 'temp_align.txt' file
		clean_up('temp_align.txt')
	
def main ():
	
	# check which method needs to be used to retrieve the representative sequence
	if args.s == 'random':
		get_rand_seq(extract_seq(args.i), size_filt(args.c, args.m), args.o)
	elif args.s == 'consensus':
		get_cons_seq(extract_seq(args.i), size_filt(args.c, args.m), False, args.o)
	elif args.s == 'combined':
		get_cons_seq(extract_seq(args.i), size_filt(args.c, args.m), args.r, args.o)	

if __name__ == "__main__":
    main()


