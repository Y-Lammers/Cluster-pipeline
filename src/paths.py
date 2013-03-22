#!/usr/bin/env python

# The paths.py script searches for the filepaths of the programs
# that will be used by this pipeline, the filepaths will be saved to
# the paths.txt file and can be manually edited if needed.

import sys

def search(program):
	# import modules to run bash commands and read the shell information
	from subprocess import call
	import os

	# get all files related to the program
	file_list = []
	for (paths, dirs, files) in os.walk('/'):
	    for file in files:
		if program == file.split('/')[-1][:len(program)]:
        	    file_list.append(os.path.join(paths, file))

	file_list.sort(key = len)
	
	# run each program in order to determine the correct path
	for item in file_list:
		try:
			p = call([item])#, stdout=open(os.devnull, 'w'))
			return item
		except:
			pass

	return 'no path found'
	
def main ():
	
	# check if there is already a paths.txt file present, if yes skip the search
	try:
		path_file = open(sys.argv[1] + 'paths.txt', 'r')
	except:
		path_file = open(sys.argv[1] + 'paths.txt', 'w')
		
		# look for the program and save the results
		for path in ['usearch', 'cd-hit', 'uclust', 'makeblastdb', 'blastn', 'tgicl', 'octu.pl']:
			filepath = search(path)
			if path == 'octu.pl': path = 'octupus'
			path_file.write(path + '\t' + filepath + '\n')

if __name__ == "__main__":
    main()
