# The pipeline.py script takes the user input and controls the various subscripts
# and programs such as blast+ and muscle.


# import the argparse module to handle the input commands
import argparse
import sys

# get the commandline arguments for the various input file, settings and output files
parser = argparse.ArgumentParser(description = 'Pipeline to process 454 reads, clusters and identifies')

parser.add_argument('--input_file', metavar='fasta files', type=str, 
			help='Enter the 454 sequence fasta file(s)', nargs='+')
parser.add_argument('--out_dir', metavar='output directory', type=str, 
			help='Enter the output directory (full path)')
parser.add_argument('--pipeline', metavar='pipeline', type=str, 
			help='The way the 454 reads will be processed (only relevant if multiple input files are used): combine (input files are combined in a single file) / cluster (input files are tagged and clustered together, bassed on tags the origin of reads in clusters can be traced) / test (run the cluster step and write the cluster info to a output file) (default: combine)', default='combine')
parser.add_argument('--trim', metavar='trim input sequences', type=str, 
			help='Trim the input sequences yes / no (default: no)', default='no')
parser.add_argument('--trim_ref', metavar='reference based trimming', type=str, 
			help='Enter the reference file for the sequence trimming (by default no reference is used)', default='no')
parser.add_argument('--save_no_match', metavar='save unligned seqs', type=str, 
			help='Save the sequences that cannot be aligned yes / no (default: yes)', default='yes')
parser.add_argument('--program', metavar='cluster program', type=str, 
			help='The cluster program that will be used to cluster the 454 reads: uclust / cdhit (default: uclust)', default='uclust')
parser.add_argument('--similarity', metavar='sequene similarity for clustering', type=str, 
			help='Sequence similarity threshold used for clustering (default: 0.97)', default='0.97')
parser.add_argument('--blast', metavar='blast method', type=str, 
			help='The blast method used for identifying the reads: genbank / Local (default: genbank', default='genbank')
parser.add_argument('--reference', metavar='reference files', type=str, 
			help='Reference file(s) used for the local blast search', nargs='+')
parser.add_argument('--blast_header', metavar='blast fasta headers', type=str, 
			help='does the reference fasta file contain blast headers (indicated by |\'s in the name) yes / no (default" no)', default='no')			
parser.add_argument('--pick_rep', metavar='otu sequence picking', type=str, 
			help='Method how the OTU representative sequence will be picked: random / consensus / combined (default: random)', default='random')
parser.add_argument('--min_size', metavar='minimum OTU size', type=int, 
			help='minimum size for an OTU to be analyzed (default: 10)', default=10)
parser.add_argument('--rand_cons', metavar='# random sequences used for the consensus', type=int, 
			help='The number of random sequences that will be pick from an OTU to make a consensus sequence (default: 10)', default=10)
args = parser.parse_args()

def get_path (pipeline):
	from subprocess import Popen, PIPE
	
	# get full path to pipeline dir
	path = Popen(['pwd'], stdout=PIPE)
	path = path.communicate()[0].replace('\n','')

	if '/home/' in path: path += '/' + '/'.join(pipeline.split('/')[:-1])
	else: path = '/' + '/'.join(pipeline.split('/')[:-1])

	return path

def check_dir (out_dir):
	import os

	# check if dir exists, if not make dir
	try:
		os.stat(out_dir)
	except:
		os.mkdir(out_dir)

def tag (pipe_path, fasta_files, out_dir):
	from subprocess import Popen, PIPE	

	# tag the input files
	tag = Popen(['python', (pipe_path + 'tag_fasta_files.py'), '-i'] + fasta_files + ['-o', out_dir], stdout=PIPE)
	tag_files = tag.communicate()[0].split('\n')[:-1]

	return tag_files

def combine (fasta_files, file_name, out_dir):
	from subprocess import call
	
	# catenate multiple fasta files into a single file for clustering
	command = 'cat ' + ' '.join(fasta_files) + ' > ' + out_dir + file_name
	#p = call(['cat'] + fasta_files + [('>' + out_dir + 'combined_fasta_file.fasta')])
	p = call(command, shell=True)

	return (out_dir + file_name)
	
def trim (pipe_path, fasta_file, reference_file, save, out_dir):
	from subprocess import call
	
	out_file = out_dir + 'trimmed_sequences.fasta'
	
	# run the trim_sequence.py script
	p = call(['python', (pipe_path + 'trim_sequence.py'), '-i', fasta_file, '-o', out_file, '-r', reference_file, '-s', save])

	return out_file

def cluster (fasta_file, similarity, program, out_dir):
	from subprocess import call
	
	# run the QIIME pick_otus.py script for fasta sequences and settings
	p = call(['pick_otus.py', '-i', fasta_file, '-o', out_dir, '-m', program, '-s', similarity])
	
	return
	
def cluster_stat (pipe_path, cluster_file, clust_time, out_dir):
	from subprocess import call
	
	# set the outputfile
	out_file = out_dir + 'cluster_stats.csv'
	
	# used the cluster_stat.py script to retrieve cluster information for the
	# cluster file
	p = call(['python', (pipe_path + 'cluster_stat.py'), '-c', cluster_file, '-t', clust_time, '-o', out_file])
	
def pick_rep_seq (pipe_path, fasta_file, cluster_file, method, min_size, rand, out_dir):
	from subprocess import call

	# get rep sequence from cluster file
	output_file = out_dir + 'clust_rep.fasta'
	if method == 'random': proc = call(['python', (pipe_path + 'pick_otu_rep.py'), '-i', fasta_file, '-o', output_file, '-c', cluster_file, '-m', str(min_size)])
 	if method == 'consensus': proc = call(['python', (pipe_path + 'pick_otu_rep.py'), '-i', fasta_file, '-o', output_file, '-c', cluster_file, '-m', str(min_size), '-s', method])
	if method == 'combined': proc = call(['python', (pipe_path + 'pick_otu_rep.py'), '-i', fasta_file, '-o', output_file, '-c', cluster_file, '-m', str(min_size), '-s', method, '-r', str(rand)])

	return output_file

def genbank_blast (pipe_path, fasta_file, out_dir):
	from subprocess import call
	
	# blast the sequences against the genbank database with the blast.py script
	output_file = out_dir + 'blast_result.csv'
	p = call(['python', (pipe_path + 'blast.py'), '-i', fasta_file, '-o', output_file])
		
	return output_file

def local_blast (pipe_path, fasta_file, reference, blast_header, out_dir):
	from subprocess import call
	
	# blast the sequences against a local blast database with the custom_blast_db.py script
	output_file = out_dir + 'blast_result.csv'
	p = call(['python', (pipe_path + 'custom_blast_db.py'), '-i', fasta_file, '-o', output_file, '-r', reference, '-b', blast_header])

	return output_file

def compare_cluster (pipe_path, cluster_file, blast_file, tag_file, min_size, out_dir):
	from subprocess import call

	# when multiple fasta files where clusted, see where the sequences in a cluster come from
	# this is done with the cluster_freq.py script
	output_file = out_dir + 'cluster_input_seqs.csv'
	p = call(['python', (pipe_path + 'cluster_freq.py'), '--cluster_file', cluster_file, '--output_file', output_file,
			'--tag_file', tag_file, '--blast', blast_file, '--min_size', str(min_size)])	

	return output_file

def main ():
	# import time module to bench the cluster run
	import time
	
	# get pipeline path
	pipe_path = get_path(sys.argv[0])
	
	# check / make the output directory
	if args.out_dir[-1] != '/': out_dir = args.out_dir + '/'
	else: out_dir = args.out_dir
	check_dir(out_dir)
		
	# check if there are multiple files that might need tagging
	if len(args.input_file) > 1:
		print('Tagging input files')
		taged_files = tag(pipe_path, args.input_file, out_dir)
		print('Merging tagged files')
		fasta_file = combine(taged_files, 'combined_fasta_file.fasta', out_dir)
	else:
		fasta_file = args.input_file[0]
	
	# trim sequences if option is selected
	if args.trim != 'no':
		print('Trimming sequences')
		fasta_file = trim(pipe_path, fasta_file, args.trim_ref, args.save_no_match, out_dir)
	
	# cluster the fasta file with the desired settings
	print('Clustering sequence file')
	time1 = time.time()
	cluster(fasta_file, args.similarity, args.program, out_dir)	
	cluster_file = out_dir + '.'.join(fasta_file.split('.')[:-1]).split('/')[-1] + '_otus.txt'
	
	# get the cluster information
	cluster_stat(pipe_path, cluster_file, time.strftime('%H:%M:%S', time.gmtime(int(time.time() - time1))), out_dir)
	
	# continue with the analysis or stop if the pipeline was only used for testing
	if args.pipeline == 'test': return
	
	# pick representative sequence for each cluster
	print('Picking representative sequences for clusters')
	rep_seq = pick_rep_seq(pipe_path, fasta_file, cluster_file, args.pick_rep, args.min_size, args.rand_cons, out_dir)
	
	# identify the clusters
	print('Identifying clusters')
	if args.blast == 'genbank':
		iden_file = genbank_blast(pipe_path, rep_seq, out_dir)
	else: 
		# check if there are multiple reference files
		if len(args.reference) > 1:
			print('Merging reference files')
			reference = combine(args.reference, 'combined_reference_file.txt', out_dir)
		else:
			reference = args.reference[0]
		
		iden_file = local_blast(pipe_path, rep_seq, reference, args.blast_header, out_dir)
	
	# combine cluster and identification files (only when multiple fasta files are used
	# and the --pipeline parameter is set to cluster
	if args.pipeline == 'cluster': compare_cluster(pipe_path, cluster_file, iden_file, (out_dir + 'tag_file.txt'), args.min_size, out_dir)
	
if __name__ == "__main__":
    main()
    
