from Bio import SeqIO
import sys

print(sys.argv)
if len(sys.argv) > 1:
    infile = sys.argv[1]
else:
    infile = "./old/MA-sequences-2-toy.fasta"

if len(sys.argv)> 2:
    size = int(sys.argv[2])
else:
    size = 250

records = list( SeqIO.parse(infile, "fasta") )

for i in range(0, len(records), size):
    outfile = infile[:-6]+'-{:03d}.fasta'.format(i//size)
    SeqIO.write(records[i:min(len(records), i+size)], outfile, "fasta")
