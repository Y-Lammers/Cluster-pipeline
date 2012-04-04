# parse the cluster results, so that for each cluster the sequences / genus and species can be listed

import argparse

# get arguments
parser = argparse.ArgumentParser(description = 'Combine the cluster and sequence files')

parser.add_argument('--cluster_file', metavar='cluster file', type=str, 
			help='enter the cluster file')
parser.add_argument('--output_file', metavar='output file', type=str, 
			help='the output file path')
parser.add_argument('--min_size', metavar='minimum OTU size', type=int, 
			help='minimum size for an OTU to be analyzed (default: 10)', default=10)
args = parser.parse_args()

def add_dic (dictionary, value):
	
	# count value
	try:
		dictionary[value] += 1
	except:
		dictionary[value] = 1

	return dictionary

def sorted_list (dictionary):
	
	# return the items in the dictionary, sorted by abundance (descending)
	item_list = ['\t'.join([k, str(v)]) for [k, v] in dictionary.items()]
	item_list.sort(reverse=True)
	
	return item_list
	

def write_output (output, output_file):
	
	#open the file
	writefile = open(output_file, 'a')
	writefile.write('\t'.join(output))
	writefile.close()

def prepare_output (cluster_info, output_file):
	
	genus, species, sequence = {}, {}, []
	
	for item in cluster_info:
		sequence.append(item[1])
		genus = add_dic(genus, item[2][1])
		species = add_dic(species, ' '.join(item[2][1:3]))
		
	genus_list = sorted_list(genus)
	species_list = sorted_list(species)
	clust = 'bla'
	
	while len(sequence) != 0:
		if len(cluster_info) > 0 and clust == 'bla':
			clust = cluster_info.pop()[0]
		else:
			clust = ''
		if len(sequence) > 0:
			seq = sequence.pop()
		else:
			seq = ''
		if len(genus_list) > 0:
			gen = genus_list.pop()
		else:
			gen = '\t'
		if len(species_list) > 0:
			spec = species_list.pop()
		else:
			spec = '\t'
		print('\n'.join([clust, seq, gen, spec]))
		try:
			write_output([clust, seq, gen, spec, '\n'], output_file)
		except:
			pass
			#write_output([cluster_info.pop()[0], sequence.pop(), genus_list.pop(), species_list.pop(), '\n'], output_file)

def extract_cluster (cluster_file, min_size, output_file):
	
	# get a list that contains the sequences present in each cluster
	cluster_list, header = [line.replace('\n','').split('\t') for line in open(cluster_file, 'r')], 'yes'
	
	# parse through the cluster list
	for cluster in cluster_list:
		if len(cluster[1:]) >= min_size:
			seqs = [[cluster[0],seq,seq.split('_')[-3:]] for seq in cluster[1:]]
			if header == 'yes':
				write_output(['cluster', 'fasta file', 'genus', '# genus', 'species', '# species', '\n'], output_file)
				header = 'no'
			prepare_output(seqs, output_file)
			
def main ():
	
	# get the cluster information
	extract_cluster(args.cluster_file, args.min_size, args.output_file)
	
if __name__ == "__main__":
    main()	
							
