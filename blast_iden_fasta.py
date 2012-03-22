# add the blast hits to the fasta header of the sequences

# import the argparse module to handle the input commands
import argparse

parser = argparse.ArgumentParser(description = 'Append the blast hits to the fasta header')

parser.add_argument('-i', metavar='fasta file', type=str, 
			help='Enter the representative fasta file')
parser.add_argument('-b', metavar='blast file', type=str, 
			help='The blast.csv produced by the pipeline')
parser.add_argument('-t', metavar='taxonomic information', type=str, 
			help='Is there taxonomic information included in the blast.csv file \'yes / no\'(yes when the genbank or unite database was used, default: no)', default='no')
parser.add_argument('-o', metavar='output fasta file', type=str, 
			help='Enter the output fasta file')
args = parser.parse_args()

def extract_seq (seq_file):
	# The SeqIO module is imported from Biopython to ease the parsing of fasta files
	from Bio import SeqIO
	
	seq_dic = {}

	# Create a dictionary of fasta sequences with the headers as keys
	for seq_record in SeqIO.parse(seq_file, 'fasta'):
		seq_dic[seq_record.id] = seq_record.seq
	
	# return the dictonary containing the sequences from the list of input files
	return seq_dic
	
def extract_blast (blast_file, tax):
	# retrieve the identifications from the blast file
	
	lines = [line.split('\t') for line in open(blast_file, 'r') if 'Percentage matched' not in line]
	blast_dic = {}
	
	# parse through the blast hits
	for line in lines:
		# add addtional taxonomic information if specified
		if tax == 'yes':
			blast_dic[line[0]] = line[12:15]
		else:
			blast_dic[line[0]] = [line[1]]
	
	return blast_dic
	
def write_results (seq, out_path):
	# import a set of Biopython modules for sequence writing
	from Bio import SeqIO

	# write a fasta sequence to the output file
	out_file = open(out_path, 'a')
	SeqIO.write(seq, out_file, 'fasta')
	out_file.close()	

def replace_header (seq_dic, blast_dic, tax, out_path):
	# import a set of Biopython modules to edit the headers
	from Bio import SeqIO
	from Bio.SeqRecord import SeqRecord
	import random, string
	
	# parse through the sequence dictionary
	for seq in seq_dic:
		header, tag = '' , ''.join(random.choice(string.ascii_uppercase + string.digits) 
				for x in range(6)) # generate the actual tag
		# check if an identification is available in the blast_dic
		if seq in blast_dic:
			try:
				header = 'cluster#:' + seq.split('_cluster_#:')[1].split('_')[0] + '_' + '_'.join(blast_dic[seq])			
			except:
				header = tag + '_' + '_'.join(blast_dic[seq])
		else:
			try:
				header = 'cluster#:' + seq.split('_cluster_#:')[1].split('_')[0] + '_no_identification_found'
			except:				
				header = tag + '_' + seq + '_no_identification_found'
		seq_rec = SeqRecord(seq_dic[seq], id=header, description='')
		write_results(seq_rec, out_path)

def main ():
	
	# get the fasta and blast information
	sequence_dic = extract_seq(args.i)
	blast_dic = extract_blast(args.b, args.t)
	
	# add the identifications to the blast headers and write the results
	replace_header(sequence_dic, blast_dic, args.t, args.o)
	
if __name__ == "__main__":
    main()


		
