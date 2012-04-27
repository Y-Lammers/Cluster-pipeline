# The custom_blast_db.py script takes a fasta file and tries to identify
# each read by blasting it against a reference set of sequences. The script
# uses the blast+ tool to setup a reference database (makeblastdb) and then
# identify the fasta sequences with the blastn tool.


# import the argparse module to handle the input commands

import argparse, sys

# get the arguments that provide the input fasta and reference file, together with the output file.
parser = argparse.ArgumentParser(description = 'create a custom blast database and run a fasta file against it')

parser.add_argument('-i', metavar='fasta file(s)', type=str, 
			help='enter the fasta file that needs to be identified')
parser.add_argument('-o', metavar='output file', type=str, 
			help='enter the output file')
parser.add_argument('-r', metavar='reference file', type=str, 
			help='enter the reference file on wich the blast database will be based')
parser.add_argument('-u', metavar='UNITE db', type=str, 
			help='UNITE databased used as reference (\'yes\' includes the taxonomic and accession numbers) yes / no (default" no)', default='no')			
parser.add_argument('-a', metavar='accession codes', type=str,
			help='Does the fasta reference contain accession codes in the header (indicated by \'|\'s) yes / no (default no)', default='no')
args = parser.parse_args()

def get_path ():
	# get the filepaths to the blast+ programs
	path_file = open('/'.join(sys.argv[0].split('/')[:-1])+'/paths.txt', 'r')
	paths, path_dic = [line.split('\t') for line in path_file], {}
	
	for line in paths:
		path_dic[line[0]] = line[1].replace('\n','')
	
	return path_dic

def make_db (ref_fasta_file, ncbi_path):
	# import module that allows the makeblastdb program to be run
	from subprocess import call
	
	#make a reference database with the makeblastdb tool
	p = call([ncbi_path, '-in', ref_fasta_file, 
			'-dbtype', 'nucl', '-parse_seqids', '-hash_index', '-max_file_sz', '5GB'])

def run_db (sequence_fasta_file, ncbi_path, db_name, outfile):
	# import module that allows the blastn tool to be run
	from subprocess import call

	#run the fasta file against the reference database with the blastn tool
	print('Blasting query file against the database . . .')
	p = call([ncbi_path, '-query', sequence_fasta_file, '-db',
			db_name, '-out', outfile, '-evalue', '10', '-outfmt', '10',
			'-max_target_seqs', '1'])
	
def clean_up (db_name, u, a):
	# import module that can run bash commands
	from subprocess import call

	#clean up database files with the 'rm' command
	if u != 'no' or a != 'no':
		cmd = 'rm ' + db_name + '*'	
	else:
		cmd = 'rm ' + db_name + '.*'
	p = call(cmd, shell = True)

def parse_blast_result (csv_path, UNITE, fasta_path):
	# parse the blast results to add headers and normalize the output so it is
	# similar to the blast.py output.

	# import modules for fasta sequence handling
	#from Bio import SeqIO

	lines = [line for line in open(csv_path, 'r')]
	#seq_list = [seq for seq in SeqIO.parse(fasta_path, 'fasta')]
	
	# write the header
	csvfile = open(csv_path, 'w')
	if UNITE != 'no':
		csvfile.write('\t'.join(['Query','Sequence','Percentage matched','length match','mismatches','gaps','query start','query end','subject start','subject end','e-value','bitscore','accession', 'genus', 'species', 'taxonomy\n']))
	else:
		csvfile.write('\t'.join(['Query','Sequence','Percentage matched','length match','mismatches','gaps','query start','query end','subject start','subject end','e-value','bitscore\n']))
	
	# add quotes around the blast hit, to simplify the importation of the csv file into spreadsheat programs
	for line in lines:
		line = line.split(',')
		line[-10] = line[-10].replace('.', ',')
		info = '\t'.join(line[-11:]).replace('\n', '')
		blast = ','.join(line[:-(len(line)-1)])
		if UNITE != 'no':
			header = line[-11].split('_')
			try:
				if 'uncultured' in header[1]:
					tax = '\t'.join([header[0], header[1].split('-')[1], header[1].split('-')[0], header[-1].replace('-', ' ')])
				else:
					tax = '\t'.join([header[0], header[1].split('-')[0], header[1].split('-')[1], header[-1].replace('-', ' ')])
			except:
				tax = '\t'.join(['','','',''])
			csvfile.write(blast + '\t' + info + '\t' + tax + '\n')
		else:
			csvfile.write(blast + '\t' + info + '\n')
		# if line[0] in seq_list: seq_list.del(line[0])
	
	# for item in seq_list:
	#	csvfile.write(item + '\tno identification found\nr')
		
def scan_fasta_file (fasta_path):
	# import modules for fasta sequence handling
	from Bio import SeqIO
	from Bio.SeqRecord import SeqRecord
	import random, string

	# read the fasta sequences
	seq_list = [seq for seq in SeqIO.parse(fasta_path, 'fasta')]
	out_path = '.'.join(fasta_path.split('.')[:-1]) + '_normalized.fasta'
	
	out_file = open(out_path, 'w')

	# replace the '|' signs that cause blastmakedb to crash, overwrite the old file
	for seq in seq_list:
		header = seq.description.split('|')
		try:
			while '' or '-' in header:
				if '' in header: header.remove('')
				if '-' in header: header.remove('-')
			if ' ' in header[1]: header[1] = header[1].replace(' ','-')
			else: header[1] += '-species'
			if len(header) == 2: header.append('no-taxonomy')
			else: header[2] = header[2].replace(' ','-')
		except:
			pass
		SeqIO.write(SeqRecord(seq.seq, id='_'.join(header), description=''), out_file, 'fasta')

	out_file.close()
	
	return out_path

def main ():
	# make the blast database
	ref_path = args.r

	if args.u != 'no' or args.a != 'no':
		ref_path = scan_fasta_file(args.r)
	
	path_dic = get_path()
	make_db(ref_path, path_dic['makeblastdb'])

	
	# blast the fasta file against the newly created database
	run_db(args.i, path_dic['blastn'], ref_path, args.o)

	# remove database files
	clean_up(ref_path, args.u, args.a)
	
	# process output
	parse_blast_result(args.o, args.u, args.i)	

		
if __name__ == "__main__":
    main()


