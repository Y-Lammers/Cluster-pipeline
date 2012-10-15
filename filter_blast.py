# Parse through the .csv file with blast results and filter this list based on 
# the percentage matched or the length of the blast match.

import argparse

# get arguments
parser = argparse.ArgumentParser(description = 'Filter the blast .csv file')

parser.add_argument('-b', metavar='blast .csv file', type=str, 
			help='enter the blast .csv file')
parser.add_argument('-o', metavar='output file', type=str, 
			help='enter the output file')
parser.add_argument('-p', metavar='minimum percentage match', type=float, 
			help='enter the minimum percentage match (default 0.0)', default = 0.0)
parser.add_argument('-l', metavar='minimum length', type=int, 
			help='enter the minimum length match (default 0)', default = 0)
args = parser.parse_args()

def filter_blast (blast_path, out_path, percentage, length):
	# open the blast file and the new blast.csv file
	blast_file = open(blast_path, 'r')
	blast_out = open(out_path, 'a')
	
	for line in blast_file:
		if 'Percentage matched' in line:
			blast_out.write(line)
		else:
			# check if the blasthit has the correct percentage matched
			# and length match, if yes save to the new output file
			line2 = line.split('\t')
			try:
				if float(line2[2]) >= percentage and int(line2[3]) >= length:
					blast_out.write(line)
			except:
				pass

	
	# close the files
	blast_file.close()
	blast_out.close()

def main ():
	# run the filtering script
	filter_blast(args.b, args.o, args.p, args.l)
	
if __name__ == "__main__":
    main()
	
