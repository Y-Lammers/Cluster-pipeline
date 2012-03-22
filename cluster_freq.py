# List for each cluster how many sequences there are present from each input file

import argparse

# get arguments
parser = argparse.ArgumentParser(description = 'List the fasta sequences present in each cluster')

parser.add_argument('--cluster_file', metavar='cluster file', type=str, 
			help='enter the cluster file')
parser.add_argument('--tag_file', metavar='tag file', type=str, 
			help='enter the tag file')
parser.add_argument('--output_file', metavar='output file', type=str, 
			help='the output file path')
parser.add_argument('--blast', metavar='blast file', type=str,
			help='the file with the blast identifications')
parser.add_argument('--min_size', metavar='minimum OTU size', type=int, 
			help='minimum size for an OTU to be analyzed (default: 10)', default=10)
args = parser.parse_args()

def get_tag (tag_file):
	
	# get the tags that are associated with the fasta files
	return [line.replace('\n', '').split('\t')[0::2] for line in open(tag_file, 'r')]

def get_blast (blast_file):
	
	# parse the blast file and retrieve the identications
	blast_dic = {}	
	for line in open(blast_file, 'r'):
		if 'Percentage matched' not in line:
			line = line.replace('\"','').replace('\n','').split('\t')
			for i in range(0, len(line)):
				if '_cluster_' in line[i]:
					blast_dic[line[i].split('_cluster_')[0]] = line[i:(i+4)]+line[12:16]

	return blast_dic
	
def write_result (tag_dic, tag_list, cluster, header, blast, out_path):

	# write the results
	out_file = open(out_path, 'a')
	if header == 'yes': out_file.write('Cluster\tBlast idenification\tPercentage matched\tlength match\taccession\tgenus\tspecies\ttaxonomy\t' + '\t'.join([item[0] for item in tag_list]) + '\n')
	temp = '\t'.join(blast)
	for item in tag_list:
		try: 
			temp += ('\t' + str(tag_dic[item[1]]))
		except:
			temp += ('\t0')
	out_file.write(temp + '\n')


def otu_freq_dist (otufile, tag_list, blast_dic, min_size, out_path):
	
	# retrieve the otu clusters that meet the size restriction
	header = 'yes'
	blast_key = blast_dic.keys()
	for line in open(otufile, 'r'):
		array = line.replace('\n', '').split('\t')[1:]
		if len(array) >= min_size:
			tag_dic, cluster, blast = {}, array[0], ['seq', 'No identification', '','','','','','']
			for seq in array:
				try:
					tag_dic[seq.split('_')[0]] += 1
				except:
					tag_dic[seq.split('_')[0]] = 1
				if seq in blast_key:
					blast = blast_dic[seq]
			if blast[0] == 'seq': blast[0] = seq
			write_result(tag_dic, tag_list, cluster, header, blast, out_path)
			header = 'no'

def main ():
	# run the comparison between input files
	otu_freq_dist(args.cluster_file, get_tag(args.tag_file), get_blast(args.blast), args.min_size, args.output_file)

	
if __name__ == "__main__":
    main()

