#!/usr/bin/env python

# The pipeline.py script takes the user input and controls the various subscripts
# and programs such as blast+ and muscle.


# import the argparse module to handle the input commands
import argparse, sys

# get the commandline arguments for the various input file, settings and output files
parser = argparse.ArgumentParser(description = 'Pipeline to process 454 reads, clusters and identifies')

parser.add_argument('-i', '--input_file', metavar='fasta files', dest='i', type=str, 
			help='Enter the 454 sequence fasta file(s)', nargs='+')
parser.add_argument('-o', '--output_dir', metavar='output directory', dest='o', type=str, 
			help='Enter the output directory (full path)')
parser.add_argument('-p', '--pipeline', metavar='merge / cluster (default) / test', dest='p', type=str, 
			help='How the input files will be processed, cluster (default) will cluster the files seperately, merge will combine them and test will only cluster without identifying', default='cluster')
parser.add_argument('--min_length', metavar='minimum length', type=int,
			help='Filter out sequences smaller then the minimum length', default=0)
parser.add_argument('--max_length', metavar='maximum length', type=int,
			help='Filter out sequences larger then the maximum length', default=0)
parser.add_argument('-c', '--program', metavar='usearch / tgicl / octupus (default) / cdhit', dest='c', type=str, 
			help='The cluster program that will be used to cluster the 454 reads: uclust / cdhit / usearch / usearch_old / tgicl / octupus (default: octupus)', default='octupus')
parser.add_argument('-s', '--similarity', metavar='clustering similarity', dest='s', type=str, 
			help='Sequence similarity threshold used for clustering (default: 0.97)', default='0.97')
parser.add_argument('-b', '--blast', metavar='genbank / locale', dest='b', type=str, 
			help='The blast method used for identifying the reads: genbank / Local (default: genbank', default='genbank')
parser.add_argument('-r', '--reference', metavar='reference files', dest='r', type=str, 
			help='Reference file(s) used for the local blast search', nargs='+')
parser.add_argument('--accession', action='store_true',
			help='Does the reference file contain accession codes in the fasta header (indicated by the |\'s)')			
parser.add_argument('-ba', '--BLAST_algorithm', metavar='algorithm', dest='ba', type=str, 
			help='Enter the algorithm BLAST wil use (default=blastn)', default='blastn')
parser.add_argument('-bd', '--BLAST_database', metavar='database', dest='bd', type=str,
			help = 'Enter the database BLAST wil use (default=nt)', default = 'nt')
parser.add_argument('-mb', '--megablast', dest='mb', action='store_true', 
			help = 'Use megablast, can only be used in combination with blastn')
parser.add_argument('-mi', '--min_identity', metavar='minimum identity', dest='mi', type=str, 
			help = 'Enter the minimal identity for BLAST results (default 97)', default = '97')
parser.add_argument('-mc', '--min_coverage', metavar='minimum coverage', dest='mc', type=str, 
			help = 'Enter the minimal coverage for BLAST results (default 100)', default = '100')
parser.add_argument('-me', '--max_evalue', metavar='maximum e-value', dest='me', type=str, 
			help = 'Enter the minimal E-value for BLAST results (default = 0.05)', default = '0.05')
parser.add_argument('-pr', '--pick_rep', metavar='random / consensus (defaul)', dest='pr', type=str, 
			help='Method how the OTU representative sequence will be picked *note* consensus mode is not supported by cd-hit and tgicl clustering, these settings will automatically select the random mode', default='consensus')
parser.add_argument('-ms', '--min_size', metavar='minimum OTU size', dest='ms', type=str, 
			help='minimum size for an OTU to be analyzed (default: 2)', default='2')
parser.add_argument('--cores', metavar='# cpu cores for tgicl', type=int,
			help='number of processors used for the tgicl cluster analysis (default: 1)', default=1)
parser.add_argument('-m', '--method', metavar='blast / species', dest='m', type=str, 
			help='The method used for merging the blast results (default = blast)', default='blast')			
args = parser.parse_args()


def get_path (pipeline):
	import os
	
	# get full path to pipeline dir
	path, file = os.path.split(os.path.realpath(__file__))
	
	return path + '/'


def check_dir (out_dir):
	import os

	# check if dir exists, if not make dir
	try:
		os.stat(out_dir)
	except:
		os.mkdir(out_dir)


def get_program_path (pipe_path):
	from subprocess import call
	
	path = call([(pipe_path + 'paths.py'), pipe_path])


def filter_seq (pipe_path, fasta_files):
	from subprocess import Popen, PIPE	

	# tag the input files
	filt = Popen([(pipe_path + 'filter.py'), '-i'] + fasta_files + ['-l', str(args.min_length), '-h', str(args.max_length)], stdout=PIPE)
	filter_files = filt.communicate()[0].split('\n')[:-1]

	return filter_files


def combine (fasta_files, file_name, out_dir):
	from subprocess import call
	
	# catenate multiple fasta files into a single file for clustering
	command = 'cat ' + ' '.join(fasta_files) + ' > ' + out_dir + file_name
	p = call(command, shell=True)

	return (out_dir + file_name)

def cluster (pipe_path, fasta_file, out_dir):
	from subprocess import call
	
	# set the output file
	out_file = out_dir + 'clustered'

	# run the cluster.py script for fasta sequences and settings
	p = call([(pipe_path + 'cluster.py'), '-i', fasta_file, '-o', out_file, '-p', args.c, '-s', args.s, '-c', str(args.cores)])

	cluster_file = out_dir + '.'.join(fasta_file.split('/')[-1].split('.')[:-1]) + '_otus.txt'

	# change the output format to a general QIIME-like format
	p = call([(pipe_path + 'cluster_to_txt.py'), '-c', out_file, '-o', cluster_file, '-p', args.c])

	return cluster_file
	

def cluster_stat (pipe_path, cluster_file, clust_time, out_dir):
	from subprocess import call
	import os	

	dir, file = os.path.split(cluster_file)
	out_file = os.path.join(dir, '_'.join(file.split('_')[:-1]) + '_stat.csv')

	# used the cluster_stat.py script to retrieve cluster information for the
	# cluster file
	p = call(['python', (pipe_path + 'cluster_stat.py'), '-c', cluster_file, '-t', clust_time, '-o', out_file])

	
def pick_rep_seq (pipe_path, fasta_file, cluster_file, out_dir):
	from subprocess import call
	import os

	# get rep sequence from cluster file
	dir, file = os.path.split(fasta_file)
	output_file = out_dir + file.split('.')[0] + '_rep.fasta'

	if args.pr == 'random' or args.c == 'cdhit' or args.c == 'tgicl': 
		proc = call(['python', (pipe_path + 'pick_otu_rep.py'), '-i', fasta_file, '-o', output_file, '-c', cluster_file, '-m', 'random', '-s', args.ms])
 	elif args.pr == 'consensus' and args.c != 'cdhit' and args.c != 'tgicl':
		fasta_file = '/'.join(cluster_file.split('/')[:-1]) + '/'
		if args.c == 'usearch' or args.c == 'usearch_old': fasta_file += 'clustered_cons'
		if args.c == 'octupus': fasta_file += 'octulist'
		proc = call([(pipe_path + 'pick_otu_rep.py'), '-i', fasta_file, '-o', output_file, '-c', cluster_file, '-m', 'consensus', '-s', args.ms, '-p', args.c])

	return output_file


def genbank_blast (pipe_path, fasta_file, out_dir):
	from subprocess import call
	
	# blast the sequences against the genbank database with the blast.py script
	output_file = out_dir + 'blast_result.csv'
	if args.mb == True:
		p = call([(pipe_path + 'blast.py'), '-i', fasta_file, '-o', output_file, '-ba', args.ba, '-bd', args.bd, '-mb', '-mi', args.mi, '-mc', args.mc, '-me', args.me])
	else:
		p = call([(pipe_path + 'blast.py'), '-i', fasta_file, '-o', output_file, '-ba', args.ba, '-bd', args.bd, '-mi', args.mi, '-mc', args.mc, '-me', args.me])
	return output_file


def local_blast (pipe_path, fasta_file, reference, out_dir):
	from subprocess import call
	
	# blast the sequences against a local blast database with the custom_blast_db.py script
	output_file = out_dir + 'blast_result.csv'
	if args.accession == True:	
		p = call([(pipe_path + 'local_blast.py'), '-i', fasta_file, '-o', output_file, '-r', reference, '-a', '-mi', args.mi, '-mc', args.mc, '-me', args.me])
	else:
		p = call([(pipe_path + 'local_blast.py'), '-i', fasta_file, '-o', output_file, '-r', reference, '-mi', args.mi, '-mc', args.mc, '-me', args.me])
	return output_file


def compare_runs (pipe_path, blast_files, output_dir):
	from subprocess import call

	# when multiple input files are clustered and identified the
	# compare_runs.py script summarizes the results in the
	# combined_blast.csv file
	output_file = output_dir + 'combined_blast.csv'
	
	p = call(['python', (pipe_path + 'compare_runs.py'), '-b'] + blast_files + ['-o', output_file, '-m', args.m])	

	return output_file


def main ():
	# import time module to bench the cluster run
	import time
	
	# get pipeline path
	pipe_path = get_path(sys.argv[0])
	
	# check / make the output directory
	if args.o[-1] != '/': output_dir = args.o + '/'
	else: output_dir = args.o
	check_dir(output_dir)
	
	# check if the program paths are set
	get_program_path(pipe_path)

	input_files = args.i
	
	# check if the inputfiles need filtering
	if args.min_length != 0 or args.max_length != 0:
		temp = filter_seq(pipe_path, input_files)
		input_files = temp
		
	# check if there are multiple files that might need tagging
	if args.p == 'merge':
		print('Merging tagged files')
		input_files = [combine(input_files, 'combined_fasta_file.fasta', output_dir)]
	
	blast_files = []

	# walk through the list of fasta files (only one when merged)
	# and cluster the sequences, identify the results, et cetera
	for fasta_file in input_files:
		
		# cluster the fasta file with the desired settings		
		print('Clustering sequence file %s' % fasta_file)
		
		time1 = time.time() # used to track the time needed for clustering
		cluster_file = cluster(pipe_path, fasta_file, output_dir)	
		
		# get the cluster information
		cluster_stat(pipe_path, cluster_file, time.strftime('%H:%M:%S', time.gmtime(int(time.time() - time1))), output_dir)
		
		# continue with the analysis or stop if the pipeline was only used for testing
		if args.p == 'test': continue
		
		# pick representative sequence for each cluster
		print('Picking representative sequences for clusters')
		rep_seq = pick_rep_seq(pipe_path, fasta_file, cluster_file, output_dir)
		
		# identify the clusters
		print('Identifying clusters')
		if args.b == 'genbank':
			iden_file = genbank_blast(pipe_path, rep_seq, output_dir)
		else: 
			# check if there are multiple reference files
			if len(args.reference) > 1:
				print('Merging reference files')
				reference = combine(args.reference, 'combined_reference_file.txt', output_dir)
			else:
				reference = args.reference[0]
			iden_file = local_blast(pipe_path, rep_seq, reference, output_dir)

			
		blast_files.append(iden_file)
				
	# combine cluster and identification files (only when multiple fasta files are used
	# and the --pipeline parameter is set to cluster
	if args.p == 'cluster': 
		print('Check the origin of the sequences in the clusters')
		compare_runs(pipe_path, blast_files, output_dir)
	
if __name__ == "__main__":
	main()
    
