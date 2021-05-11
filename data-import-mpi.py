import sys
import traceback
import multiprocessing

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

    print(frequency)
    mut_rate = frequency
    print(mut_rate)

    return mut_rate


def main():

    # Read Covid19 reference sequence
    ref = SeqIO.read("./data/ref_sequence.gb", "genbank")
    print('Reference Covid sequence')
    print(ref.id)
    print(repr(ref.seq))
    print(len(ref.seq))


    # Read downloaded sequence file from NCBI GenBank Virus site
    
    if len(sys.argv) > 2: 
        infile = sys.argv[1]
    else:
        infile = "./data/MA-sequences-2-toy.fasta"
    records = list( SeqIO.parse(infile, "fasta") )
    
    metadata = []
    mutarray = []
    
    import time
    start = time.time()
    for idx, record in enumerate(records):
        print('\n{} of {} records'.format(idx+1, len(records)))
        try:
            meta = get_meta_fasta(record)
            mut = mutation_array(ref.seq, record.seq)
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

if __name__ == "__main__":
    main()
