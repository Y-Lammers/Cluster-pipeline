#!/usr/bin/env python

# The cluster.py script takes the input sequence file and clusters settings
# and runs the cluster program of choice

# import the argparse module to handle the input commands
import argparse, sys

# get the 4 commandline arguments for the input files and output directory
parser = argparse.ArgumentParser(description = 'Cluster sequences')

parser.add_argument('-i', metavar='sequence file', type=str, 
			help='Enter the sequence file path')
parser.add_argument('-o', metavar='output file/dir', type=str,
			help='The output file/dir name')			
parser.add_argument('-p', metavar='program', type=str,
			help='The program used for clustering (octupus, usearch, usearch_old, cdhit, tgicl) default: octupus', default='usearch')
parser.add_argument('-s', metavar='similarity', type=str,
			help='The sequence similarity used for clustering', default='0.97')
parser.add_argument('-c', metavar='cpu cores', type=int,
			help='the number of cpu cores used for tgicl clustering (default = 1)', default=1)
args = parser.parse_args()

def get_path ():
	# get the filepaths to the cluster programs
	path_file = open('/'.join(sys.argv[0].split('/')[:-1])+'/paths.txt', 'r')
	paths, path_dic = [line.split('\t') for line in path_file], {}
	
	for line in paths:
		path_dic[line[0]] = line[1].replace('\n','')
	
	return path_dic

def command (cmd):
	# import module that can run bash commands
	from subprocess import call

	#clean up database files with the 'rm' command
	p = call(cmd, shell = True)
	
def run_usearch (usearch):
	# import module that allows the usearch tool to be run
	from subprocess import call

	p = call([usearch, '-cluster_fast', args.i, '-id', args.s, '-uc', args.o, '-consout', args.o + '_cons'])

def run_usearch_old (usearch):
	# import module that allows the usearch tool to be run
	from subprocess import call
	from Bio import SeqIO
	from Bio.SeqRecord import SeqRecord
	from Bio.Seq import Seq

	# remove N characters that can interfere with the -consout flag in usearch

	seq_dic = SeqIO.index(args.i, "fasta")
	out_file = open(args.i + '.sorted', 'w')
	for item in seq_dic:
		temp = seq_dic[item]
		temp = SeqRecord(Seq(str(temp.seq).replace('N','')), id=temp.id, description='')
		SeqIO.write(temp, out_file, "fasta")
	out_file.close()
	
	p = call([usearch, '-sort', args.i + '.sorted', '-output', args.i + '.sorted'])
	p = call([usearch, '-cluster', args.i + '.sorted', '-uc', args.o, '-id', args.s, '-consout', args.o + '_cons'])
	command('rm ' + args.i + '.sorted')

def run_tgicl (tgicl, output_dir):
	# import module that allows the tgicl tool to be run
	from subprocess import call

	# run the tgicl program
	p = call([tgicl, '-F', args.i, '-p', args.s, '-c', str(args.c)])

	# set the correct output directory (by default the desired files remain in the directory
	# where the tgicl cluster program is runned
	if output_dir[-1] != '/': output_dir = '/'.join(output_dir.split('/')[:-1]) + '/'


	run_directory = '/'.join(sys.argv[0].split('/')[:-1]) + '/'

	# a list of command that will move and remove all cluster files created by tgicl
	cmd_list = ['mv ' + run_directory + '*_cl_clusters ' + output_dir + 'clustered_cl_clusters',
			'mv ' + run_directory + '*.singletons ' + output_dir + 'clustered.singletons',
			'rm ' + sequence_file + '.*', 'rm ' + run_directory + 'masked.lst',
			'rm ' + run_directory + '*tgicl*', 'rm ' + run_directory + '*.log', 'rm ' + run_directory + '*_cl_*',
			'rm -rf ' + run_directory + 'asm_*']

	for cmd in cmd_list:
		command(cmd)

def run_octupus (octu, similarity, output_dir):
	# import module that allows the octupus tool to be run
	from subprocess import call
	
	run_directory = '/'.join(sys.argv[0].split('/')[:-1]) + '/'

	# convert the similarity scores to the format used by octupus
	if similarity == '1.0': similarity = '1'
	if '0.' in similarity: similarity.replace('0.','')

	# set the output directory
	if output_dir[-1] != '/': output_dir = '/'.join(output_dir.split('/')[:-1]) + '/'

	# run the octupus cluster program
	p = call([octu, args.i, similarity, '0', output_dir])

	# rename the octupus cluster file so it can be used by the other scripts downstream
	cmd_list = ['mv \"' + output_dir + 'octuall.seq\" \"' + output_dir + 'clustered\"', 
			'rm ' + run_directory + '*.log']
	for cmd in cmd_list:
		p = call(cmd, shell = True)


def run_cdhit (cdhit):
	# import module that allows the cdhit tool to be run
	from subprocess import call
	
	#p = call([cdhit, '-i', ('\"' + sequence_file + '\"'), '-o',  output_file, '-c', similarity, '-d', '0'])
	p = call([cdhit, '-i', args.i, '-o',  args.o, '-c', args.s, '-d', '0'])
	
def main ():
	
	# get the cluster program paths
	paths = get_path()
	
	if args.p == 'usearch':
		run_usearch(paths['usearch'])
	if args.p == 'usearch_old':
		run_usearch_old(paths['usearch'])
	elif args.p == 'cdhit':
		run_cdhit(paths['cd-hit'])
	elif args.p == 'tgicl':
		run_tgicl(paths['tgicl'], args.o)
	elif args.p == 'octupus':
		run_octupus(paths['octupus'], args.s, args.o)
		
if __name__ == "__main__":
    main()

