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
args = parser.parse_args()

def writeresult (result, out_path, header):
	# write the blast results to the output file

	outfile = open(out_path, 'a')
	# checks if there is already a header present, if not the header will be writen
	# to the output file
	if header == 'yes':
		outfile.write('Blast hit\tSequence\tPercentage matched\tlength match\tmismatches\tgaps\tquery start\tquery end\tsubject start\tsubject end\te-value\tbitscore\tspecies\ttaxonomy\n')
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
	
	

def parse_seq_file (seq_path, out_path):
	# import the biopython module to deal with fasta parsing
	from Bio import SeqIO
	
	# a header will be written before the blast results
	header = 'yes'

	# parse the fasta file
	for seq in SeqIO.parse(seq_path, 'fasta'):
	
		previous, blastrun, count = 0, '', 0
	
		# Try to get the blast results
		while blastrun == '':
			blastrun = blast_sequence(seq)
			count += 1
			# stop trying after a 100 timeouts
			if count % 100 == 0:
				blastrun = 'error'
				print('timeout sequence %s' % seq.id)
		# if a blasthit could be retrieved from the sever
		if blastrun != 'error' and blastrun.alignments != []:
			for alignment in blastrun.alignments:
				for hsp in alignment.hsps:

					# select the best blast hit			
					if previous == 0:
				
						# get the blast information that is
						# needed to write to the results
						match_length = len(hsp.match)
						mismatch = len(hsp.match.replace('|',''))
						percent_match = str(100 - ((float(mismatch) / match_length)*100))
	
						if list(hsp.frame)[0] == 1: query_end = str(hsp.query_start + match_length)
						else: query_end = str(hsp.query_start - match_length)
	
						if list(hsp.frame)[1] == 1: sbjct_end = str(hsp.sbjct_start + match_length)
						else: sbjct_end = str(hsp.sbjct_start - match_length)
	
						previous = 1
	
						# get the taxonomic and species information
						tax_org = obtain_tax(alignment.title.split('|')[1])
						taxonomy, organism = tax_org[0], tax_org[1]
					
						# prepare and write the output
						output = '\t'.join([('\"' + alignment.title + '\"'), seq.id, percent_match, str(match_length), str(mismatch), str(hsp.gaps), str(hsp.query_start), query_end, str(hsp.sbjct_start), sbjct_end, str(hsp.expect), str(hsp.bits), organism, taxonomy, '\n'])
						writeresult(output, out_path, header)
						header = 'no'

		# write the 'no blast hit found' message to the output file if no
		# blast result could be obtained for a fasta sequence
		else: 
			writeresult('\t'.join(['no blast hit found', seq.id, '\n']), out_path, header)
	
parse_seq_file(args.i, args.o)

