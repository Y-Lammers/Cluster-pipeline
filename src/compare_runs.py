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
	import os

	blast_dic, blast_list = {}, []

	for blast_file in args.b:
		
		# grap the blast file name without path
		dir, temp_file = os.path.split(blast_file)
		blast_list.append(temp_file) # append to the file

		# parse through the file and obtain the blast info
		for line in open(blast_file, 'r'):
			if 'query,hit' not in line:
				blast_info = split_blast_line(line)			
				if args.m == 'blast':
					blast_dic = append_dic(blast_dic, blast_info[1], blast_info, temp_file)
				else:
					try:
						blast_dic = append_dic(blast_dic, blast_info[7], blast_info, temp_file)
					except:
						blast_dic = append_dic(blast_dic, blast_info[1], blast_info, temp_file)

	# return both the blast information and the blast files
	return [blast_dic, blast_list]

def append_dic(blast_dic, dic_key, blast_info, blast_file):
	
	# append the blast_info to the blast_dic, the blast dic
	# is a nested dic, and for each item there is a sub
	# dic with the info for each blast file

	# the number of sequences in the cluster for the blasted sequence
	seq_number = int(blast_info[0].split('length_cluster:')[1][:-1])

	# temp_info contains the blast hit + species and taxonomy information
	try:
		temp_info = ['\"' + blast_info[1] + '\"'] + blast_info[7:9]
	except:
		temp_info = ['\"' + blast_info[1] + '\"']

	# if the blast key is not present in the dic: add it 
	if dic_key not in blast_dic:
		blast_dic[dic_key] = {blast_file: [temp_info, seq_number]}
	else:	
		# try to add the seq number, if not possible: create a new record
		try:
			blast_dic[dic_key][blast_file][1] += seq_number
		except:
			blast_dic[dic_key][blast_file] = [temp_info, seq_number]		

	# return the updated dictionary	
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


def parse_results (blast_results):

	# For each blast hit the general info is obtained (blast hit + species)
	# and the number of sequences in the clusters that match this hit
	blast_dic, blast_list = blast_results[0], blast_results[1]	
	
	for hit in blast_dic:
		temp_result = []
		for item in blast_dic[hit]:
			temp_result += blast_dic[hit][item][0]
			break
		for item in blast_list:
			if item in blast_dic[hit]:
				temp_result.append(str(blast_dic[hit][item][1]))
			else:
				temp_result.append('0')

		# append the results for the hit to the output file
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

