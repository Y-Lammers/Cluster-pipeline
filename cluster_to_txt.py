# This script converts the output files cluster programs such as uclust
# cd hit to a generalized output format that is also compatible with the
# QIIME pipeline.

# import the argparse module to handle the input commands
import argparse

# get the 3 commandline arguments for the input files and output directory
parser = argparse.ArgumentParser(description = 'convert cluster file')

parser.add_argument('-c', metavar='cluster file', type=str, 
			help='enter the cluster file path')
parser.add_argument('-o', metavar='output file', type=str,
			help='The output file name')			
parser.add_argument('-p', metavar='program used', type=str,
			help='The program used for clustering (uclust/usearch/cdhit')
args = parser.parse_args()

def get_uc_cluster (cluster_file):
	# parse through the uc file and get the cluster information
	cluster_dic = {}

	for line in open(cluster_file, 'r'):
		line = line.replace('\n','').split('\t')
		# get the seeds and sequences who match these seeds
		if line[0] == 'S':
			cluster_dic[int(line[1])+1] = [line[8]]
		elif line[0] == 'H':
			cluster_dic[int(line[1])+1] += [line[8]]
	
	return cluster_dic

def get_cd_cluster (cluster_file):
	# parse through the cd hit file and get the cluster information
	
	cluster_dic, cluster = {}, 0
	
	for line in open(cluster_file + '.clstr', 'r'):
		line = line.replace('\n','')
		# get the cluster number if a new cluster starts
		if line[0] == '>':
			cluster = int(line.split(' ')[1])+1
			cluster_dic[cluster] = []
		# expand the cluster with the sequences in it
		else:
			seq = line.split('>')[1].split('...')[0]
			cluster_dic[cluster] += [seq]
				
	return cluster_dic

def get_octupus_cluster (cluster_file):
	# parse the octupus cluster file and get the cluster information
	# cluster file = octuall.seq

	cluster_dic, cluster = {}, 0

	for line in open(cluster_file, 'r'):
		line = line.replace('\n','')
		if '*' in line:
			cluster = int(line.replace('*octu',''))
		if '>' in line:
			try:
				cluster_dic[cluster].append(line[1:])
			except:
				cluster_dic[cluster] = [line[1:]]
	return cluster_dic			

def get_tgicl_cluster (cluster_file):
	# parse the tgicl cluster file and get the cluster information
	
	cluster_dic, cluster = {}, 0

	# parse the (non-singleton) cluster file and add the clusters to the dictionary
	for line in open(cluster_file + '_cl_clusters', 'r'):
		line = line.replace('\n','').split('\t')
		if '>' in line[0]:
			cluster = int(line[0][3:])
		else:
			cluster_dic[cluster] = line

	cluster_singleton = cluster + 1

	# parse the singleton cluster file and add the singletons to the dictionary
	for line in open(cluster_file + '.singletons', 'r'):
		cluster_dic[cluster] = [line.replace('\n','')]
		cluster += 1
	
	return cluster_dic	

def write_output (cluster_dic, output_file):
	# convert the clusters to the desired .txt output and
	# write them to the output file
	
	# open the output file
	output = open(output_file, 'w')
	
	for i in range(1, (len(cluster_dic)+1)):
		cluster =  '\t'.join([str(i)] + cluster_dic[i]) + '\n'
		output.write(cluster)
	
	output.close()
		
def main ():
	
	if args.p == 'uclust' or args.p == 'usearch':
		# get the uc file and save it in the .txt format
		write_output(get_uc_cluster(args.c), args.o)
	
	elif args.p == 'cdhit':
		# get the clstr file and save it in the .txt format
		write_output(get_cd_cluster(args.c), args.o)

	elif args.p == 'octupus':
		# get the octuall.seq file and save it in the .txt format
		write_output(get_octupus_cluster(args.c), args.o)

	elif args.p == 'tgicl':
		# get the tgicl cluster and singelton files, extract the clusters and merge the results
		write_output(get_tgicl_cluster(args.c), args.o)
		
if __name__ == "__main__":
    main()

		
