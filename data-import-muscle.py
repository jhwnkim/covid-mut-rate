from Bio import SeqIO

import traceback

# Read Covid19 reference sequence
ref = SeqIO.read("./data/ref_sequence.gb", "genbank")
print('Reference Covid sequence')
print(ref.id)
print(repr(ref.seq))
print(len(ref.seq))

def get_meta_fasta(record):
    import re

    tokens = re.split(r'\|', record.description)

    metadata = {
        'id': record.id,
        'collect-date': tokens[-1],
        'country': tokens[-2]
    }

    print(metadata)
    return metadata


from Bio import Align
from Bio.Align.substitution_matrices import Array

from Bio.Align.Applications import MuscleCommandline

import subprocess
# muscle_bin = r"/home/gridsan/jhwnkim/git-repo/tools/muscle3.8.31/src/muscle"
# muscle_in = r"/home/gridsan/jhwnkim/tmp/muscle_in.fasta"
# muscle_out= r"/home/gridsan/jhwnkim/tmp/muscle_out.fasta"
muscle_bin = "./muscle3.8.31_i86win32.exe"
muscle_in = "./tmp/muscle_in.fasta"
muscle_out= "./tmp/muscle_out.fasta"

def mutation_array(seq1, seq2): # pass as SeqRecord
    SeqIO.write([seq1, seq2], muscle_in, "fasta")
    #muscle_cline = MuscleCommandline(muscle_bin, input=muscle_in)
    #print(muscle_cline)
    #child = subprocess.Popen(muscle_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    muscle_cline = [muscle_bin, "-in", muscle_in]
    child = subprocess.Popen(muscle_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    print(muscle_cline)
    # aligned_seq = SeqIO.read(muscle_out, "fasta")
    aligned_seq = list(SeqIO.parse(child.stdout, "fasta"))

    se1 = aligned_seq[0].seq
    se2 = aligned_seq[1].seq

    frequency = Array("ACGTN-", dims=2)
    for c1, c2 in zip(se1, se2):
        frequency[c1, c2] += 1

    print(frequency)
    mut_rate = frequency
    #print(mut_rate)

    return mut_rate


# Read downloaded sequence file from NCBI GenBank Virus site

infile = "./data/old/MA-sequences-2-toy.fasta"
records = list( SeqIO.parse(infile, "fasta") )

metadata = []
mutarray = []

import time
start = time.time()
for idx, record in enumerate(records):
    print('\n{} of {} records'.format(idx+1, len(records)))
    try:
        meta = get_meta_fasta(record)
        mut = mutation_array(ref, record)
    except:
        print(traceback.format_exc())
    else:
        metadata.append(meta)
        mutarray.append(mut)



dates = []
mutarray_avg = []
ids = []

for idx, rec in enumerate(metadata):
    if len(mutarray_avg) ==0 or rec['collect-date'] > dates[-1]:
        dates.append(rec['collect-date'])
        ids.append([rec['id']])
        mutarray_avg.append(mutarray[idx])

    else:
        for i in range(len(mutarray_avg)):
            if rec['collect-date']< dates[i]:
                dates.insert(i, rec['collect-date'])
                ids.insert(i,[rec['id']])
                mutarray_avg.insert(i, mutarray[idx])
                break
            elif rec['collect-date'] == dates[i]:
                ids[i].append(rec['id'])
                mutarray_avg[i] += mutarray[idx]
                break

# Divide mutation rate by counts and convert to float Array
mutvec_out = []
for idx, idlist in enumerate(ids):
    mutarray_avg[idx] = mutarray_avg[idx]/len(idlist)

    print(dates[idx])
    print(idlist)
    print(mutarray_avg[idx])

    mutvec =  [ \
		mutarray_avg[idx]["A", "C"], \
		mutarray_avg[idx]["A", "T"], \
		mutarray_avg[idx]["A", "G"], \
		mutarray_avg[idx]["C", "A"], \
		mutarray_avg[idx]["C", "T"], \
		mutarray_avg[idx]["C", "G"], \
		mutarray_avg[idx]["T", "A"], \
		mutarray_avg[idx]["T", "C"], \
		mutarray_avg[idx]["T", "G"], \
		mutarray_avg[idx]["G", "A"], \
		mutarray_avg[idx]["G", "C"], \
		mutarray_avg[idx]["G", "T"]]
    mutvec_out.append([float(x) for x in mutvec])

# print('{:e}'.format(mutarray_avg[0]['C', 'T']))
# Save to file
import pandas as pd

# outfile = './data/MA-sequences-1-toy1.csv'
# outfile = './data/MA-sequences-2-toy.csv'
outfile = infile[:-5] + '.csv'

df = pd.DataFrame({"Dates": dates})
df = pd.concat( [df, \
	pd.DataFrame(mutvec_out, columns=['A->C', 'A->T', 'A->G', \
									'C->A', 'C->T', 'C->G', \
									'T->A', 'T->C', 'T->G', \
									'G->A', 'G->C', 'G->T'])], axis=1)
df = pd.concat( [df, \
	pd.DataFrame({"N": [len(idlist) for idlist in ids]})], axis=1)
print(df)
df.to_csv(outfile, index=False)

print('Run time took {}'.format(time.strftime("%H:%M:%S", time.gmtime(time.time()-start))))

'''
To read data
# load file with pandas

import pandas as pd

outfile = './data/MA-sequences-1-toy1.csv'

df = pd.read_csv(outfile)

# convert to list and numpy array
dates = df['Dates'].values.tolist() # in strings
mutrates = df.iloc[:,1:13].to_numpy()


print(df)
'''
