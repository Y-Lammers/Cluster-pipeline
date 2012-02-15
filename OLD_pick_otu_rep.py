import sys, random

def size_filt(otufile, minsize):
	
	count, arraydic = 1, {}	

	for line in open(otufile, 'r'):
		array = line.replace('\n', '').split('\t')[1:]
		if len(array) >= minsize: 
			arraydic[count] = array
		count += 1

	return arraydic

def extract_seq(seqfile):

	seqdic, header, seq = {}, '', ''

	for line in open(seqfile, 'r'):
		if line[0] == '>':
			if seq != '':
				seqdic[header] = seq
				seq = ''
			header = line.split(' ')[0][1:]
		else:
			seq += line.replace('\n', '')
	seqdic[header] = seq

	return seqdic
	
def get_rep_seq(seqdic, arraydic, max_result, outpath):
	from random import shuffle

	outfile = open(outpath, 'w')

	keys = arraydic.keys()
	
	shuffle(keys)

	for item in keys[:max_result]:
		header = random.choice(arraydic[item])
		outfile.write('>' + header + '\tcluster #: ' + str(item) + '\tlength cluster: ' + str(len(arraydic[item])) + '\n' + seqdic[header] + '\n')

get_rep_seq(extract_seq(sys.argv[1]), size_filt(sys.argv[2], int(sys.argv[3])), int(sys.argv[4]), sys.argv[5])

