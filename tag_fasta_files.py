# The tag_fasta_files.py script provides a fasta file with 'tags' for each
# sequence header, based on these tags the input files from which sequences
# came from can be traced after clustering.


# import the argparse module to handle the input commands
import argparse

# get the 2 commandline arguments for the input files and output directory
parser = argparse.ArgumentParser(description = 'Add unique tags to fasta files')

parser.add_argument('-i', metavar='fasta files', type=str, 
			help='enter the fasta file(s)', nargs='+')
parser.add_argument('-o', metavar='output dir', type=str,
			help='the directory in which the files will be placed (default: same folder as input files)', default='')			
args = parser.parse_args()


def retrieve_fasta_files (fasta_file_list, out_dir):
	# import the Biopython modules to both parse and write fasta files
	from Bio import SeqIO
	from Bio.SeqRecord import SeqRecord
	# import the random and string modules to generate a random 'tag' to append to the sequences
	import random, string
	

	tag_list, tag = [], ''

	# go through the list of input fasta files
	for fasta_path in fasta_file_list:

		# for each fasta file a 6 character unique tag is generated
		while tag in tag_list or tag == '': # check if the tag is unique
			tag = ''.join(random.choice(string.ascii_uppercase + string.digits) 
				for x in range(6)) # generate the actual tag
		else:
			# store the tag in a list to prevent the same tag being used twice
			tag_list.append(tag)

		# define outputfiles and output directory
		if out_dir != '':
			# if no output directory is specified, the tagged files will be placed
			# in the same directory as the input files
			output_path = (out_dir + '.'.join(fasta_path.split('/')[-1].split('.')[:-1]) + 
					'_tag_' + tag + '.fasta')
		else:
			# if a output directory is specified the files will be stored here
			output_path = ('.'.join(fasta_path.split('.')[:-1]) + '_tag_' + 
					tag + '.fasta')
		print(output_path) # the output path is printed, the paths will be
				   # used later in the pipeline
		
		# The Biopython modules is used to parse through the fasta file
		for seq in SeqIO.parse(fasta_path, 'fasta'):
			# for each sequence the tag is added to the sequence header (seq.id)
			new_seq = SeqRecord(seq.seq, id=(tag + '_' + seq.id), description='')
			
			# the resulting sequence is written to the output file
			output_file = open(output_path, 'a')
			SeqIO.write(new_seq, output_file, 'fasta')
			output_file.close()
		
		# write file names and tags to the tag_file.txt file based on this info the 
		# origin of sequences can be traced to their original fasta files
		if out_dir != '': tag_file = open((out_dir + 'tag_file.txt'), 'a')
		else: tag_file = open('tag_file.txt', 'a')
		tag_file.write(fasta_path + '\t' + output_path + '\t' + tag + '\n')
		tag_file.close()

retrieve_fasta_files(args.i, args.o)



