# The cluster.py script takes the input sequence file and clusters settings
# and runs the cluster program of choice

# import the argparse module to handle the input commands
import argparse, sys

# get the 4 commandline arguments for the input files and output directory
parser = argparse.ArgumentParser(description = 'Cluster sequences')

parser.add_argument('-i', metavar='sequence file', type=str, 
			help='Enter the sequence file path')
parser.add_argument('-o', metavar='output file', type=str,
			help='The output file name')			
parser.add_argument('-p', metavar='program', type=str,
			help='The program used for clustering (uclust, usearch, cdhit) default: usearch', default='usearch')
parser.add_argument('-s', metavar='similarity', type=str,
			help='The sequence similarity used for clustering', default='0.97')
args = parser.parse_args()

def get_path ():
	# get the filepaths to the cluster programs
	path_file = open('/'.join(sys.argv[0].split('/')[:-1])+'/paths.txt', 'r')
	paths, path_dic = [line.split('\t') for line in path_file], {}
	
	for line in paths:
		path_dic[line[0]] = line[1].replace('\n','')
	
	return path_dic

def clean_up (file_name):
	# import module that can run bash commands
	from subprocess import call

	#clean up database files with the 'rm' command
	cmd = 'rm ' + file_name
	p = call(cmd, shell = True)
	
def run_uclust (uclust, sequence_file, similarity, output_file):
	# import module that allows the uclust tool to be run
	from subprocess import call
	p = call([uclust, '--sort', sequence_file, '--output',  sequence_file + '.sorted'])
	p = call([uclust, '--input', sequence_file + '.sorted', '--uc', output_file, '--id', similarity])
	clean_up(sequence_file + '.sorted')

def run_usearch (usearch, sequence_file, similarity, output_file):
	# import module that allows the usearch tool to be run
	from subprocess import call
	
	p = call([usearch, '-sort', sequence_file, '-output',  sequence_file + '.sorted'])
	p = call([usearch, '-cluster', sequence_file + '.sorted', '-uc', output_file, '-id', similarity])
	clean_up(sequence_file + '.sorted')

def run_cdhit (cdhit, sequence_file, similarity, output_file):
	# import module that allows the cdhit tool to be run
	from subprocess import call
	
	p = call([cdhit, '-i', sequence_file, '-o',  output_file, '-c', similarity, '-d', '0'])
	
def main ():
	
	# get the cluster program paths
	paths = get_path()
	
	if args.p == 'usearch':
		run_usearch(paths['usearch*'], args.i, args.s, args.o)
	elif args.p == 'cdhit':
		run_cdhit(paths['cd-hit'], args.i, args.s, args.o)
	elif args.p == 'uclust':
		run_uclust(paths['uclust'], args.i, args.s, args.o)
		
if __name__ == "__main__":
    main()

