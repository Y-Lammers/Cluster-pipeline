import sys

'''from Bio import SeqIO

# create a dictionary of fasta sequences with the headers as keys
seq_dic = SeqIO.index(sys.argv[1], 'fasta')

#for item in seq_dic: print(item)

keys = seq_dic.keys()

for item in keys: print(seq_dic[item].id)

'''
from subprocess import Popen, PIPE	
#muscle directory
muscle_path = Popen(['find', '/home' ,'-name', 'muscle*'], stdout=PIPE)
muscle_path = muscle_path.communicate()[0].replace('\n', '')
from subprocess import call
p = call([muscle_path, '-in', sys.argv[1], '-out', 'temp_align.txt', '-quiet'])
from Bio import AlignIO
alignment = AlignIO.read('temp_align.txt', 'fasta')
#print(alignment)
#gap_consensus(self, threshold=0.7, ambiguous='X',
#consensus_alpha=None, require_multiple=0
from Bio.Align import AlignInfo
#summary_align = AlignInfo.SummaryInfo(alignment)
#consensus = summary_align.dumb_consensus()
consensus = AlignInfo.SummaryInfo(alignment).dumb_consensus()
print(consensus)


##help(AlignInfo.SummaryInfo)

#p = call(['rm', 'temp_align.txt'])

