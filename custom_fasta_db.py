import sys, subprocess

def paths():
	from subprocess import Popen, PIPE	
	
	#working directory
	work_path = Popen(['pwd'], stdout=PIPE)
	work_path = work_path.communicate()[0].replace('\n', '')

	#ncbi bin directory
	ncbi_path = Popen(['find', '/home' ,'-name', 'makeblastdb'], stdout=PIPE)
	ncbi_path = ncbi_path.communicate()[0].split('makeblastdb')[0]

	return(work_path, ncbi_path)

def make_db (ref_fasta_file, ncbi_path):
	from subprocess import call
	
	#make database
	db_name = (ref_fasta_file.split('.')[0] + '_ref_data_base')

	p = call([(ncbi_path + 'makeblastdb'), '-in', ref_fasta_file, '-title',
			db_name, '-dbtype', 'nucl', '-parse_seqids', '-input_type',
			'fasta', '-hash_index'])
	print('data base created')


def run_db (sequence_fasta_file, ncbi_path, db_name, outfile):
	from subprocess import call

	#run query
	p = call([(ncbi_path + 'blastn'), '-query', sequence_fasta_file, '-db',
			db_name, '-out', outfile, '-evalue', '10', '-outfmt', '10'])
	
	print('finished blasting the sequences')

	#clean up database files
	cmd = 'rm ' + db_name + '.*'	
	p = call(cmd, shell = True)
	print(cmd + ' removed')
	
make_db(sys.argv[1], paths()[1])
run_db(sys.argv[2], paths()[1], sys.argv[1], sys.argv[3])
print('done')


