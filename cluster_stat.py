# The cluster_stat.py script goes through the cluster output files
# and gathers information about the clusterings such as the number
# of clusters, size of the clusters and number of used sequences.
# the results will be writen to a csv output file.


# import the argparse module to handle the input commands
import argparse

# get the commandline arguments for the input files and output directory and user settings
parser = argparse.ArgumentParser(description = 'Gather information about the formed clusters')

parser.add_argument('-c', metavar='cluster file', type=str, 
			help='enter the cluster (otu) file')
parser.add_argument('-t', metavar='time', type=str, 
			help='time it took for the cluster process to complete', default='')
parser.add_argument('-o', metavar='output file', type=str, 
			help='enter the output file')

args = parser.parse_args()

def parse_otu (otufile):
	# Parse through the OTU file to get the information
	
	otu_seq_dic, sequences, otu_size_dic = {}, 0, {}

	# parse through the OTU file
	for line in open(otufile, 'r'):
		# split each OTU into a set of fasta headers that make up
		# each OTU
		array = line.replace('\n', '').split('\t')
		otu_seq_dic[array[0]] = array[1:]
		sequences += len(array[1:])
		try:
			otu_size_dic[len(array[1:])] += 1
		except:
			otu_size_dic[len(array[1:])] = 1
	
	return [otu_seq_dic, sequences, otu_size_dic]
	
def write_results (otu_info_list, otufile, time, outfile):
	# get the results and write them to the output file
	
	# open the output file
	outcsv = open(outfile, 'w')
	
	# write basic stats
	outcsv.write(otufile + '\nNumber of sequences\t' + str(otu_info_list[1]) + '\nNumber of clusters\t' + 
			str(len(otu_info_list[0])) + '\nCluster time\t' + time + '\nSize of cluster\tnumber of clusters for size\n')
	
	# for each cluster size, print the results
	for size in sorted(otu_info_list[2].iterkeys()):
		outcsv.write(str(size) + '\t' + str(otu_info_list[2][size]) + '\n')
	
	outcsv.close()
	
def main ():
	
	write_results(parse_otu(args.c), args.c, args.t, args.o)	

if __name__ == "__main__":
	main()

