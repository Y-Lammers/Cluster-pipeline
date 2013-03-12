#!/usr/bin/env python

# go throught the blast results of a n-number of cluster runs
# and combine the results in one overview table

import argparse

# get arguments
parser = argparse.ArgumentParser(description = 'List the fasta sequences present in each cluster')

parser.add_argument('-o', metavar='output file', type=str, 
			help='the output csv file')
parser.add_argument('-b', metavar='blast files', type=str,
			help='the file with the blast identifications', nargs='+')
parser.add_argument('-m', metavar='blast/species', type=str,
			help='Merger the blast files based on species or exact blast hit', default='blast')
args = parser.parse_args()


def read_blast ():
	
	# parse through the list of blast files and obtain
	# the blast information for each file. The info is stored
	# in a dictionary and the keys are either the raw blast hits
	# or the species for the blast hit
	blast_dic = {}

	for blast_file in args.b:
		for line in open(blast_file, 'r'):
			if 'query,hit' not in line:
				blast_info = split_blast_line(line)			
				if args.m == 'blast':
					blast_dic = append_dic(blast_dic, blast_info[1], blast_info, blast_file)
				else:
					blast_dic = append_dic(blast_dic, blast_info[7], blast_info, blast_file)
	
	return blast_dic

def append_dic(blast_dic, dic_key, blast_info, blast_file):
	
	# append the blast_info to the blast_dic, the blast dic
	# is a nested dic, and for each item there is a sub
	# dic with the info for each blast file

	# the number of sequences in the cluster for the blasted sequence
	seq_number = blast_info[0].split('length_cluster:')[1][:-1]

	try:
		blast_dic[dic_key][blast_file] = [['\"' + blast_info[1] + '\"'] + blast_info[7:9], seq_number]
	except:
		blast_dic[dic_key] = {blast_file: [['\"' + blast_info[1] + '\"'] + blast_info[7:9], seq_number]}

	return blast_dic
	

def split_blast_line (line):
	
	# split the blast into a list, while dealing with
	# the comma's that might appear in the blast output

	# remove the newline
	line = line.strip()

	# split the line based on the quotes around the text
	comma_split = line.split('\"')

	# split the remaining entries based on the comma, and return this list	
	return comma_split[1::2] + comma_split[4].split(',')[1:]

def blast (line):
	line = line.split('\"')
	print line[2]

def species (line):
	pass

def parse_results (blast_dic):
	
	for hit in blast_dic:
		temp_result = []
		for item in blast_dic[hit]:
			temp_result += blast_dic[hit][item][0]
			break
		for item in args.b:
			if item in blast_dic[hit]:
				temp_result.append(blast_dic[hit][item][1])
			else:
				temp_result.append('0')
		write_results(','.join(temp_result), 'a')


def write_results (result, mode):
	
	# write the results to the output file
	out_file = open(args.o, mode)
	out_file.write(result + '\n')
	out_file.close()


def main ():

	blast_results = read_blast()
	write_results('blast result,species,taxonomy,' + ','.join(args.b), 'w')
	parse_results(blast_results)


if __name__ == "__main__":
	main()

