# The blast.py script takes a fasta file and blasts each fasta file 
# against the Genbank database. The best blast results (lowest e-value)
# is saved along with bit-score, percentage match, species and the
# full taxonomic information.


# import the argparse module to handle the input commands
import argparse

# get the commandline arguments that specify the input fastafile and the output file
parser = argparse.ArgumentParser(description = ('retrieve blast and taxonomic information for a fasta file'))

parser.add_argument('-i', metavar='fasta file', type=str, 
			help='enter the fasta file')
parser.add_argument('-o', metavar='output file', type=str, 
			help='enter the output file')
parser.add_argument('-t', metavar='number of blast threads', type=int, 
			help='enter the number of blast threads, warning a large number of threads can slow down the search (max: 10)', default=4)			
args = parser.parse_args()

def writeresult (result, out_path):
	# write the blast results to the output file
	outfile = open(out_path, 'a')
	# checks if there is already a header present, if not the header will be writen
	# to the output file
	outfile.write(result)
	outfile.close()

def obtain_tax (code):
	# a module from Biopython is imported to connect to the Entrez database
	from Bio import Entrez
		
	taxonomy, organism = '', ''

	# based on the taxid the species and taxonomy are retrieved from the Entrez database
	try:
		Entrez.email = "quick@test.com"
		handle = Entrez.efetch(db="nucleotide", id= code, retmode="xml")
		record = Entrez.read(handle)
		taxonomy = record[0]["GBSeq_taxonomy"]
		organism = record[0]["GBSeq_organism"]
		handle.close()
	except:
		pass

	return(taxonomy, organism)

def blast_sequence (seq):
	# The blast modules are imported from biopython
	from Bio.Blast import NCBIWWW, NCBIXML

	# Blast the sequence against the NCBI nucleotide database
	try:
		result_handle = NCBIWWW.qblast('blastn', 'nt', seq.format('fasta'))
		blast_result = NCBIXML.read(result_handle)
		return blast_result
	except:
		return ''		
	
def get_blast_align (seq):
	# import the biopython module to deal with fasta parsing
	from Bio import SeqIO
	import time
	
	# Try to get the blast results
	blastrun, time1 = blast_sequence(seq), time.time()

	# Keep trying to get a result for 5 minutes if no result was obtained
	while blastrun == '' and time.time()-time1 < 300:
		align = blast_sequence(seq)

	# return a timeout message if no result could be obtained
	if blastrun == '':
		blastrun = 'error'
		print('timeout sequence %s' % seq.id)		

	return blastrun
	
def get_output (hsp, seq, alignment):
	import time
	
	# get the blast information that is
	# needed to write to the results
	match_length = len(hsp.match)
	mismatch = len(hsp.match.replace('|',''))
	percent_match = str(100 - ((float(mismatch) / match_length)*100))

	if list(hsp.frame)[0] == 1: query_end = str(hsp.query_start + match_length)
	else: query_end = str(hsp.query_start - match_length)

	if list(hsp.frame)[1] == 1: sbjct_end = str(hsp.sbjct_start + match_length)
	else: sbjct_end = str(hsp.sbjct_start - match_length)

	# Keep trying to get a taxonomic and species information
	# if there is no result after 5 minutes no taxonomic information will be included
	tax_org, time1 = ['',''], time.time()
	while tax_org[0] == '' and time.time()-time1 < 300:
		tax_org = obtain_tax(alignment.title.split('|')[1])

	taxonomy, organism = tax_org[0], tax_org[1]
				
	# prepare the output
	output = '\t'.join([('\"' + alignment.title + '\"'), seq.id, percent_match, 
				str(match_length), str(mismatch), str(hsp.gaps), str(hsp.query_start), 
				query_end, str(hsp.sbjct_start), sbjct_end, str(hsp.expect), 
				str(hsp.bits), organism, taxonomy, '\n'])

	return output


def parse_blast_align (seq, out_path):
	# import the biopython module to deal with fasta parsing
	from Bio import SeqIO
	
	# get the best blasthit
	previous, align = 0, get_blast_align(seq)

	if align != 'error' and align.alignments != []:
		for alignment in align.alignments:
			for hsp in alignment.hsps:

				# select the best blast hit			
				if previous == 0:
					
					# write the results to the output file
					writeresult(get_output(hsp, seq, alignment), out_path)
					previous = 1
	# write the 'no blast hit found' message to the output file if no
	# blast result could be obtained for a fasta sequence
	else: 
		writeresult('\t'.join(['no blast hit found', seq.id, '\n']), out_path)
		
	return
		
def parse_seq_file (seq_path, threads, out_path):
	# import the biopython module to deal with fasta parsing
	# and the multiprocessing module to run multiple blast threads
	from Bio import SeqIO
	import multiprocessing
	
	# parse the fasta file
	seq_list = [seq for seq in SeqIO.parse(seq_path, 'fasta')]
	
	# blast each sequence in the seq_list list
	procs = []
	while len(seq_list) > 0 or len(procs) > 0:
		# start the maximum number of threads
		while len(procs) < threads:
			try:
				p = multiprocessing.Process(target=parse_blast_align, args=(seq_list.pop(0), out_path,)) 
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
	
def parse_blast_result (csv_path):
	# parse the blast results to add headers and normalize the output
	lines = [line for line in open(csv_path, 'r')]
	
	# write the header
	csvfile = open(csv_path, 'w')
	csvfile.write('\t'.join(['Query','Sequence','Percentage matched','length match',
			'mismatches','gaps','query start','query end','subject start',
			'subject end','e-value','bitscore','species','taxonomy\n']))
	
	# add quotes around the blast hit, to simplify the importation of the csv file into spreadsheat programs
	for line in lines: 
		csvfile.write(line)

def main ():
	# check the if the number of threads is larger then 10, if yes set the
	# maximum number of threads to 10
	if args.t > 10: parse_seq_file(args.i, 10, args.o)
	else: parse_seq_file(args.i, args.t, args.o)
	
	# normalize the output
	parse_blast_result(args.o)

if __name__ == "__main__":
    main()
