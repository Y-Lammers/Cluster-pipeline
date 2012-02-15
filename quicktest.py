import sys

def cluster (fasta_file):
	from subprocess import call
	
	# run the pick_otus.py script with the selected fasta files and cluster method
	command = 'cat ' + ' '.join(fasta_file) + ' > bla.txt'
	print(command)
#	proc = call(['cat', fasta_file[0], fasta_file[1], '> test.txt'])
	proc = call(command, shell=True)

#cluster([sys.argv[1], sys.argv[2]])

a = 'abcdefghij'
'''b = {}

for i in range(0,len(a)):
	b[a[i]*2] = i

print(b)

if 'aa' in b: print('ja')
'''

for let in [let for let in a if let in 'ab']:
	print(let)

