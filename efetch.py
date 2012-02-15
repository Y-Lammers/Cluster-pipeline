from Bio import Entrez, SeqIO

code = 301040922

handle = Entrez.efetch(db="nucleotide", id= code, rettype="gb")
record = SeqIO.read(handle, "genbank")
handle.close()
print(record)
