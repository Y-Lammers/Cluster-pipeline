#!/usr/bin/env python

# The blast.py script takes a fasta file and blasts each fasta file 
# against the Genbank database. The best blast results (lowest e-value)
# is saved along with bit-score, percentage match, species and the
# full taxonomic information.


# import the argparse module to handle the input commands
import argparse

# get the commandline arguments that specify the input fastafile and the output file
parser = argparse.ArgumentParser(description = ('retrieve blast and taxonomic information for a fasta file'))

parser.add_argument('-i', '--input_file', metavar='fasta file', dest='i', type=str, 
			help='enter the fasta file')
parser.add_argument('-o', '--output_file', metavar='output file', dest='o', type=str, 
			help='enter the output file')
parser.add_argument('-ba', '--BLAST_algorithm', metavar='algorithm', dest='ba', type=str, 
			help='Enter the algorithm BLAST wil use (default=blastn)', default='blastn')
parser.add_argument('-bd', '--BLAST_database', metavar='database', dest='bd', type=str,
			help = 'Enter the database BLAST wil use (default=nt)', default = 'nt')
parser.add_argument('-mb', '--megablast', dest='mb', action='store_true', 
			help = 'Use megablast, can only be used in combination with blastn')
parser.add_argument('-mi', '--min_identity', dest='mi', type=int, 
			help = 'Enter the minimal identity for BLAST results', default = 97)
parser.add_argument('-mc', '--min_coverage', dest='mc', type=int, 
			help = 'Enter the minimal coverage for BLAST results', default = 100)
parser.add_argument('-me', '--max_evalue', dest='me', type=float, 
			help = 'Enter the minimal E-value for BLAST results', default = 0.05)
args = parser.parse_args()


def obtain_tax (code):
	# a module from Biopython is imported to connect to the Entrez database
	from Bio import Entrez
		
	taxon_info = []

	# based on the taxid the species and taxonomy are retrieved from the Entrez database
	try:
		Entrez.email = "quick@test.com"
		handle = Entrez.efetch(db="nucleotide", id= code, retmode="xml")
		record = Entrez.read(handle)
		taxon_info.append(record[0]["GBSeq_organism"]) # organism
		taxon_info.append(record[0]["GBSeq_taxonomy"]) # taxonomy
		handle.close()
	except:
		pass

	# return the taxonomy information in list form
	return taxon_info


def blast_bulk (sequences, thread_number):

	# The blast modules are imported from biopython
	from Bio.Blast import NCBIWWW, NCBIXML
	from Bio import SeqIO
	import os

	# grap the file path to the input file, here the temporary fasta files will be created
	dir, file = os.path.split(args.i)

	# create the list where all the blast results are stored in
	blast_list = []

	# create the temporary file
	temp_file_path = os.path.join(dir, (str(thread_number) + 'temp.fasta'))
	temp_file = open(temp_file_path, 'w')

	# fill the temp file with sequences
	SeqIO.write(sequences, temp_file, 'fasta')
	temp_file.close()

	# read the temp fasta file
	temp_fasta_file = open(temp_file_path, 'r')
	fasta_handle = temp_fasta_file.read()

	# blast the temporary file, and save the blasthits in the blast_list
	result_handle = NCBIWWW.qblast(args.ba, args.bd, fasta_handle, megablast=args.me, hitlist_size=1)

	blast_list += [item for item in NCBIXML.parse(result_handle)]

	# remove the temporary file		
	os.remove(temp_file_path)

	# return the filled blast hit
	return blast_list

	
def parse_blast_align (sequences, thread, mode):
	# import the biopython module to deal with fasta parsing

	blast_list = blast_bulk(sequences, thread)
	count = 1	

	# parse though the blast hits
	for blast_result in blast_list:
		for alignment in blast_result.alignments:
			for hsp in alignment.hsps:
	            		
				# calculate the %identity
		            	identity = float(hsp.identities/(len(hsp.match)*0.01))

				# grab the genbank number
				gb_num = alignment.title.split('|')[1]

				# an alignment needs to meet 3 criteria before 
				# it is an acceptable result: above the minimum 
				# identity, minimum coverage and E-value
			
				# create containing the relevant blast results
				# pass this list to the filter_hits function to
				# filter and write the blast results
				filter_hits([('\"' + blast_result.query + '\"'), ('\"' + alignment.title + '\"'), gb_num, str(identity),
						str(blast_result.query_length), str(hsp.expect), str(hsp.bits)], mode, thread, count)
				count += 1


def filter_hits (blast, mode, thread, count):
	
	# filter the blast hits, based on the minimum
	# identity, minimum coverage, e-value and the user blacklist
	if float(blast[3]) >= args.mi and int(blast[4]) >= args.mc and float(blast[5]) <= args.me:
		taxon = obtain_tax(blast[2])
		results = blast+taxon
		
		# write the results
		write_results(','.join(results), mode)


def write_results (result, mode):
	
	# write the results to the output file
	out_file = open(args.o, mode)
	out_file.write(result + '\n')
	out_file.close()

	
def parse_seq_file ():
	# import the biopython module to deal with fasta parsing
	# and the multiprocessing module to run multiple blast threads
	from Bio import SeqIO
	import multiprocessing
	import time
	
	# parse the fasta file
	seq_list, sub_list = [seq for seq in SeqIO.parse(args.i, 'fasta')], []
	
	# blast each sequence in the seq_list list
	procs, count, threads = [], 1, 10
	while len(seq_list) > 0 or len(procs) > 0:
		# start the maximum number of threads
		while len(procs) < threads and len(seq_list) > 0:
			if len(seq_list) >= 50:
				sub_list = seq_list[:50]
				seq_list = seq_list[50:]
			else:
				sub_list = seq_list
				seq_list = []
			try:
				p = multiprocessing.Process(target=parse_blast_align, args=(sub_list, count, 'a',)) 
				procs.append([p, time.time()])
				p.start()
				count+=1
			except:
				break

		# check when a thread is done, remove from the thread list and start
		# a new thread
		while len(procs) > 0:
			for p in procs:
				if p[0].is_alive() == False: 
					p[0].join()
					procs.remove(p)
				# time-out after 30 minutes
				elif time.time() - p[1] > 10800:
					p[0].terminate()
					procs.remove(p)
			break	
	

def main ():
	# create a blank result file and write the header
	#header = 'query,hit,accession,identity,hit length,e-value,bit-score,taxon id,genbank record species,CITES species,CITES info,NCBI Taxonomy name,appendix'
	header = 'query,hit,accession,identity,hit length,e-value,bit-score,taxon id,species,taxonomy'
	write_results(header, 'w')

	# Blast the sequences	
	parse_seq_file()


if __name__ == "__main__":
	main()
