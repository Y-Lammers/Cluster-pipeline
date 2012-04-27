# the trim_sequence.py script will load a set of sequences and
# trim each sequence to a certain length based on a reference sequence
# or sequences from the input file itself. Trimed sequences will
# be saved to a output fasta file, sequences that cannot be trimmed
# will be discarded or saved as well.


# import the argparse module to handle the input commands
import argparse, sys

# get the commandline arguments that specify the input fastafile and the output file
parser = argparse.ArgumentParser(description = ('Trim sequences to make these easier to cluster, trimming is based on the alignements of certain markers that the region of interst'))

parser.add_argument('-i', metavar='fasta file', type=str, 
			help='enter the fasta file')
parser.add_argument('-o', metavar='output file', type=str, 
			help='enter the output file')
parser.add_argument('-r', metavar='reference file', type=str, 
			help='enter a reference file, when no reference is specified sequences from the input fasta file will be used as reference', default='no')
parser.add_argument('-s', metavar='save untrimmed', type=str, 
			help='save untrimmed sequences yes/no (default: yes', default='yes')
parser.add_argument('-n', metavar='marker samplings', type=int, 
			help='The number of sample alignments that will be ran to determine the marker sequences (default: 1000)', default=1000)			
parser.add_argument('-t', metavar='number of threads', type=int, 
			help='The number of threads that will be used for trimming (default = 2)', default=2)			
args = parser.parse_args()

def extract_seq (seq_file):
	# The SeqIO module is imported from Biopython to ease the parsing of fasta files
	from Bio import SeqIO
	
	seq_dic = {}

	# parse through a list of fasta file
	# Create a dictionary of fasta sequences with the headers as keys
	for seq_record in SeqIO.parse(seq_file, 'fasta'):
		seq_dic[seq_record.id] = seq_record.seq
	
	# return the dictonary containing the sequences from the list of input files
	return seq_dic

def find_muscle_path ():
	# find the path to the muscle program
	path_file = open('/'.join(sys.argv[0].split('/')[:-1])+'/paths.txt', 'r')
	muscle_path = [line.split('\t')[1].replace('\n','') for line in path_file if 'muscle' in line]
		
	return muscle_path[0]

def clean_up (file_path):
	from subprocess import call
	
	# used to remove the temporary files created by the muscle alignments
	p = call(['rm', file_path])

def write_results (sequence, header, out_path):
	# import a set of Biopython modules for sequence handling and writing
	from Bio import SeqIO
	from Bio.SeqRecord import SeqRecord	

	# write a fasta sequence to the output file
	try:
		out_file = open(out_path, 'a')
		SeqIO.write(SeqRecord(sequence, id=header, description=''), out_file, 'fasta')
		out_file.close()
	except:
		pass

def run_alignment (fasta_path):
	# import the module to start other programs
	from subprocess import call
	
	# get the muscle sequence alignment path
	muscle_path = find_muscle_path()
	file_name = fasta_path.replace('.fasta', '')

	# align the newly created 'temp_in.fasta' file with muscle
	p = call([muscle_path, '-in', (file_name + '.fasta'), '-out', (file_name + '_aligned.txt'), '-quiet'])
	# remove the 'temp_in.fasta' file since this file is no longer needed
	clean_up((file_name + '.fasta'))
	
	return
	
def get_fragment (file_path):
	# Imports the necessairy biopython modules to handle multiple sequence
	# alignment files and fasta files
	from Bio import AlignIO
	from Bio.Align import AlignInfo	

	# read the sequence alignment file
	alignment, length_dic = AlignIO.read(file_path, 'fasta'), {}
	
	# remove the alignment file
	clean_up(file_path)

	# remove the gaps from the smallest sequence
	if str(alignment[0].seq).count('-') > str(alignment[1].seq).count('-'):
		splitseq = str(alignment[0].seq).split('-')
	elif str(alignment[0].seq).count('-') < str(alignment[1].seq).count('-'):
		splitseq = str(alignment[1].seq).split('-')
	else:
		splitseq = str(alignment[0].seq)

	# store the length of each framgent in the length dictionary, sort the key list to 
	# obtain the longest ungapped fragment and return this
	for fragment in splitseq:
		length_dic[len(fragment)] = fragment
	length_list = length_dic.keys()
	length_list.sort(reverse=True)

	return length_dic[length_list[0]]

def get_marker_seq (sequence_dic, reference, samples):
	# import the modules to random sample items from lists
	from random import sample
	
	# get keys for the sequence and reference dictionary, if there is no reference
	# set the sequence set will act as a reference as well.
	sequence_keys = sequence_dic.keys()
	if reference == 'no':	reference = sequence_dic
	reference_keys = reference.keys()
	
	marker_dic = {}
	
	# check if -n samples is smaller then the number of sequences, if not samples ==
	# len(sequences)
	if samples > len(sequence_keys): samples = len(sequence_keys)
	
	# sample a x number of sequences from the sequence dictionary
	for sequence in sample(sequence_keys, samples):
		# write the randomly sample sequence and reference sequence to the
		# temp_in.fasta file
		header = sample(reference_keys, 1)[0]
		file_name = '-'.join([sequence, header])
		write_results(sequence_dic[sequence], sequence, (file_name + '.fasta'))
		write_results(reference[header], header, (file_name + '.fasta'))
		
		# run the alignment for sequence file
		run_alignment(file_name + '.fasta')
		
		# obtain the largest fragment
		fragment = get_fragment((file_name + '_aligned.txt'))
		
		# get the first and last 10 bases from the fragment and store these
		# in the marker_dic, count the frequency if they appear more than once
		marker = '-'.join([fragment[:10],fragment[-10:]])
		if len(fragment) >= len(sequence)/1.33 and len(marker) == 21:
			if marker not in marker_dic:
				marker_dic[marker] = 1
			elif marker in marker_dic:
				marker_dic[marker] += 1

	# get the most abundant marker from the marker_dic and return this value
	marker_list = [[marker_dic[marker], marker] for marker in marker_dic]
	marker_list.sort(reverse=True)
	
	return (marker_list)[0][1]

def get_marker_pos (marker, file_path):
	# Imports the necessairy biopython modules to handle multiple sequence
	# alignment files and fasta files
	from Bio import AlignIO
	from Bio.Align import AlignInfo	

	# read the multiple sequence alignment file
	alignment = AlignIO.read(file_path, 'fasta')

	# remove the alignment file
	clean_up(file_path)

	return str(alignment[1].seq).find(marker)

def trim_sequence (sequence_list, markers, save, output_path):
	# import the biopython modules to deal with sequences
	from Bio.Seq import Seq
	from Bio.Alphabet import IUPAC

	marker_pos = []

	# run a alignment with the sequence for both markers to determine the trimming positions
	for marker in markers:
		if marker not in str(sequence_list[0]):
			file_name = '-'.join([sequence_list[1], marker])
			write_results(sequence_list[0], sequence_list[1], (file_name + '.fasta'))
			write_results(Seq(marker, IUPAC.unambiguous_dna), marker, (file_name + '.fasta'))

			# run the alignment for sequence file
			run_alignment(file_name + '.fasta')
			marker_pos.append(get_marker_pos(marker, (file_name + '_aligned.txt')))
		else:
			marker_pos.append(str(sequence_list[0]).find(marker))

	if marker_pos[0] < marker_pos[1] and len(sequence_list[0][marker_pos[0]:marker_pos[1]]) > 10:
		write_results(sequence_list[0][marker_pos[0]:marker_pos[1]], sequence_list[1], output_path)
	else:
		if save == 'yes' and len(sequence_list[0]) > 10:
			write_results(sequence_list[0], sequence_list[1], output_path)

	return


def trim_sequences_process (sequence_dic, marker, threads, save, output_path):
	# import the multiprocessing module to start multiple threads
	import multiprocessing
	
	sequence_list, markers, procs = [[sequence_dic[header], header] for header in sequence_dic], marker.split('-'), []

	while len(sequence_list) > 0 or len(procs) > 0:
		# start the maximum number of threads
		while len(procs) < threads and len(sequence_list) > 0:
			try:
				p = multiprocessing.Process(target=trim_sequence, args=(sequence_list.pop(0), markers, save, output_path,)) 
				procs.append(p)
				p.start()
			except:
				break
		# check when a thread is done, remove from the thread list and start
		# a new thread
		while len(procs) > 0:
			for p in procs:
				if p.is_alive() == False: 
					p.join()
					procs.remove(p)
			break	

def main ():
	# get the sequence dictionary and obtain the marker
	sequence_dic = extract_seq(args.i)

	if args.r != 'no': 
		reference = extract_seq(args.r)
	else:
		reference = 'no'
	print('obtaining marker sequences')
	marker = get_marker_seq(sequence_dic, reference, args.n)
	
	# trim the sequences
	print('trimming sequences')
	trim_sequences_process(sequence_dic, marker, args.t, args.s, args.o)
	
if __name__ == "__main__":
    main()

