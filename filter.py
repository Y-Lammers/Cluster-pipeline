# This script can parse and filter fasta input files (filter by length) or
# remove duplicate sequences

# import the argparse module to handle the input commands
import argparse

# get the 2 commandline arguments for the input files and output directory
parser = argparse.ArgumentParser(description = 'Filter input files')

parser.add_argument('-i', metavar='fasta files', type=str, 
			help='enter the fasta file(s)', nargs='+')
parser.add_argument('-o', metavar='output dir', type=str,
			help='the directory in which the files will be placed (default: same folder as input files)', default='')			
parser.add_argument('-m', metavar='min length', type=int,
			help='filter sequences based on the minimum sequence length', default=0)
parser.add_argument('-d', metavar='remove duplicate', type=str,
			help='remove duplicate sequences from the dataset (yes/no) default: no', default='no')
args = parser.parse_args()

# get sequence dictionary
def get_seq (fasta_path, length, dup):
	# import modules for fasta sequence handling
	from Bio import SeqIO
	
	seq_dic, used_dic = {}, {}
	# parse the fasta file
	for seq in SeqIO.parse(fasta_path, 'fasta'):
		# check if the sequence is longer then the threshold
		if len(seq.seq) > length:
			if dup != 'no':
				# check if the sequence isn't a duplicate
				if seq.seq not in used_dic:
					seq_dic[seq.id] = seq.seq
					used_dic[seq.seq] = 1
			else:
				seq_dic[seq.id] = seq.seq
	
	return seq_dic

def write_fasta (seq_dic, out_path):
	# import the modules for fasta sequence handling
	from Bio import SeqIO
	from Bio.SeqRecord import SeqRecord
	
	out_file = open(out_path, 'w')
	for seq in seq_dic:
		new_seq = SeqRecord(seq_dic[seq], id=(seq), description='')
		SeqIO.write(new_seq, out_file, 'fasta')
	out_file.close()
	
def main ():
	
	# filter each input file
	
	for fasta_file in args.i:
		seq_dic = get_seq(fasta_file, args.m, args.d)
		if args.o != '':
			out_path = (args.o + '.'.join(fasta_file.split('/')[-1].split('.')[:-1]) + 
					'_filtered.fasta')
		else:
			out_path = '.'.join(fasta_file.split('.')[:-1]) + '_filtered.fasta'
		# print out_path for use later in the pipeline
		print(out_path)
		write_fasta(seq_dic, out_path)
		
if __name__ == "__main__":
    main()
    	
