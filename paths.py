# The paths.py script searches for the filepaths of the programs
# that will be used by this pipeline, the filepaths will be saved to
# the paths.txt file and can be manually edited if needed.

import sys

def search(program):
	# import modules to run bash commands and read the shell information
	from subprocess import Popen, PIPE

	# use the linux 'find' command to find the program
	path = Popen(['find', '/home/', '/usr/', '/software/' ,'-name', program], stdout=PIPE)
	path = path.communicate()[0].split('\n')[0]
	
	return path
	
def main ():
	
	# check if there is already a paths.txt file present, if yes skip the search
	try:
		path_file = open(sys.argv[1] + 'paths.txt', 'r')
		pass
	except:
		path_file = open(sys.argv[1] + 'paths.txt', 'w')
		
		# look for the program and save the results
		for path in ['muscle*', 'usearch*', 'cd-hit', 'uclust', 'makeblastdb', 'blastn']:
			filepath = search(path)
			path_file.write(path + '\t' + filepath + '\n')

if __name__ == "__main__":
    main()
