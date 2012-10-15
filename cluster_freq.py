# Go through the list of clusters, and list how many sequences there are present in each cluster and
# from what input file these sequences originate. The script can merge different clusters together
# if they have the same blast hit, with this it possible to combine different markers.

import argparse

# get arguments
parser = argparse.ArgumentParser(description = 'List the fasta sequences present in each cluster')

parser.add_argument('-c', metavar='cluster txt file', type=str, 
			help='enter the cluster text file')
parser.add_argument('-f', metavar='cluster fasta file', type=str, 
			help='enter the cluster fasta file')
parser.add_argument('-t', metavar='tag file', type=str, 
			help='enter the tag file')
parser.add_argument('-o', metavar='output file', type=str, 
			help='the output file path')
parser.add_argument('-b', metavar='blast file', type=str,
			help='the file with the blast identifications')
parser.add_argument('-m', metavar='merge clusters', type=str,
			help='Merge the clusters based on the blast hits (yes/no), default = no', default='no')
parser.add_argument('-s', metavar='minimum OTU size', type=int, 
			help='minimum size for an OTU to be analyzed (default: 10)', default=10)
args = parser.parse_args()

def get_tag (tag_file):
	
	# get the tags that are associated with the fasta files
	return [line.replace('\n', '').split('\t')[0::2] for line in open(tag_file, 'r')]

def get_blast (blast_file):
	
	# parse the blast file and retrieve the identications
	blast_dic, seq_dic = {}, {}
	for line in open(blast_file, 'r'):
		if 'Percentage matched' not in line:
			line = line.replace('\"','').replace('\n','').split('\t')
			for i in range(0, len(line)):
				if '_cluster_' in line[i]:
					if len(line) > 15:
						blast_dic[line[i+1]] = line[i:(i+4)]+line[12:16]
					else:
						blast_dic[line[i+1]] = line[i:(i+4)]+['','','','']
					seq_dic[line[i].split('_cluster_')[0]] = line[i+1]

	return [blast_dic, seq_dic]
	
def get_otu (otu_file):
	
	# return a nested list with all the sequences in the OTU's
	return [line.replace('\n','').split('\t')[1:] for line in open(otu_file, 'r')]

def get_cluster_name (fasta_file):
	# import the biopython module to process the fasta file
	from Bio import SeqIO
	
	cluster_dic = {}
	
	# parse through the cluster fasta file, and retrieve the cluster name for each otu
	for seq_record in SeqIO.parse(fasta_file, 'fasta'):
		cluster_dic[seq_record.id.split('_cluster')[0]] = seq_record.id
	
	return cluster_dic

def combine (otu_list, seq_dic, cluster_dic, merge):
	
	# combine the blast and otu results, if the merge parameter is set to yes,
	# clusters will be combined if they have the same blast hit
	
	count, seq_key, combine_dic, cluster_key = 0, seq_dic.keys(), {}, cluster_dic.keys()
	
	# parse through the set otus
	for otu in otu_list:
		blast, cluster = 'No identification' + str(count), ''
		# for each sequence present in the otu
		for seq in otu:
			# check if there is a blast hit available for the sequence
			if seq in seq_key:
				blast = seq_dic[seq]
			# check if there sequence is also the name of a full otu
			if seq in cluster_key:
				cluster = [cluster_dic[seq]]
		# if the merge paramter is set to 'merge', the otu information will be
		# merged with other otus that share the same blast information
		if merge == 'yes':
			try:
				combine_dic[blast][0] += otu
				combine_dic[blast][1] += cluster
			except:
				combine_dic[blast] = [otu, cluster]
		# add a unique number to the blast hit if the data will not be merged,
		# this will prevent overwrites over other clusters that share the same
		# blast information
		else:
			combine_dic[str(count) + '_____num_____' + blast] = [otu, cluster]
		count += 1
	
	return combine_dic
			
	
def write_result (tag_dic, tag_list, header, blast, out_path):

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

def otu_freq_dist (combine_dic, tag_list, blast_dic, min_size, out_path):
	
	# retrieve the otu clusters that meet the size restriction
	header = 'yes'
	blast_key = blast_dic.keys()
	
	# go through the dictionary that contains the otu's
	for item in combine_dic:
		# check if the length of the otu is larger then the size threshold
		if len(combine_dic[item][0]) >= min_size:
			tag_dic, blast = {}, ['seq', 'No identification', '','','','','','']
			# go through the otu sequences and count based on the tag
			# how many sequences there are present from each input dataset
			for seq in combine_dic[item][0]:
				try:
					tag_dic[seq.split('_')[0]] += 1
				except:
					tag_dic[seq.split('_')[0]] = 1
			if '_____num_____' in item:
				item2 = item.split('_____num_____')[1]
			else:
				item2 = item
			# try to obtain the expanded blast description for the cluster
			if item2 in blast_key:
				blast = blast_dic[item2]
			blast[0] = '- '.join(combine_dic[item][1])
			write_result(tag_dic, tag_list, header, blast, out_path)
			header = 'no'			

def main ():
	# get the tags
	tag = get_tag(args.t)
	
	# get the blast information
	blast =  get_blast(args.b)
	
	# get otus sequences
	otu_seqs = get_otu(args.c)
	
	# get the otu names
	otu_name = get_cluster_name(args.f)
	
	# combine and export the information
	otu_freq_dist(combine(otu_seqs, blast[1], otu_name, args.m), tag, blast[0], args.s, args.o)
	
if __name__ == "__main__":
    main()

