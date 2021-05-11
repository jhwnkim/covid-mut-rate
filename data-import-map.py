import os, sys, json
import traceback

from Bio import SeqIO
from Bio import Align
from Bio.Align.substitution_matrices import Array

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


def mutation_array(seq1, seq2):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 2
    aligner.mismatch_score = -3
    aligner.open_gap_score = -7
    aligner.extend_gap_score = -2

    alignments = aligner.align(seq1, seq2)
    alignment = alignments[0]

    frequency = Array("ACGTN", dims=2)
    for (start1, end1), (start2, end2) in zip(*alignment.aligned):
        se1 = seq1[start1:end1]
        se2 = seq2[start2:end2]
        for c1, c2 in zip(se1, se2):
            frequency[c1, c2] += 1

    return frequency


# Read Covid19 reference sequence
ref = SeqIO.read("./data/ref_sequence.gb", "genbank")
print('Reference Covid sequence {}'.format(ref.id))

# Read downloaded sequence file from NCBI GenBank Virus site
sys.path.append('../')

inOutFileList = open(sys.argv[1],"r+")

import time
start = time.time()

for line in inOutFileList.readlines():

    # Get input and output file names for this iteration
    (inFile,outFile) = line.split()
    print("Reading from " + inFile + " and writing to " + outFile)

    infile = inFile
    records = list( SeqIO.parse(infile, "fasta") )

    multi_list = []

    for idx, record in enumerate(records):
        print('\n{} of {} records'.format(idx+1, len(records)))
        try:
            meta = get_meta_fasta(record)
            mut = mutation_array(ref.seq, record.seq)
        except:
            print(traceback.format_exc())
        else:
            mutvec = [ \
	    		mut["A", "C"], \
	    		mut["A", "T"], \
	    		mut["A", "G"], \
	    		mut["C", "A"], \
	    		mut["C", "T"], \
	    		mut["C", "G"], \
	    		mut["T", "A"], \
	    		mut["T", "C"], \
	    		mut["T", "G"], \
	    		mut["G", "A"], \
	    		mut["G", "C"], \
	    		mut["G", "T"]]

            # row: Date, id, country, mutation rates
            multi_list.append( [meta['collect-date'], meta['id'], meta['country']] + [float(x) for x in mutvec])

            print(multi_list[-1])
            print('{} Elapsed'.format(time.strftime("%H:%M:%S", time.gmtime(time.time()-start))))


	# sort data according to collected data starting from oldest to new
    sorted_multi_list = sorted(multi_list, key=lambda x: x[0])

    # Save to file
    import pandas as pd

    # outfile = './data/MA-sequences-1-toy1.csv'
    # outfile = './data/MA-sequences-2-toy.csv'
    outfile = infile[:-5] + 'csv'

    df = pd.DataFrame(sorted_multi_list, columns=['Dates', 'ID', 'Geolocation', \
										'A->C', 'A->T', 'A->G', \
    									'C->A', 'C->T', 'C->G', \
    									'T->A', 'T->C', 'T->G', \
    									'G->A', 'G->C', 'G->T'])
    print(df)
    df.to_csv(outfile, index=False)

    print('Run time took {}'.format(time.strftime("%H:%M:%S", time.gmtime(time.time()-start))))
