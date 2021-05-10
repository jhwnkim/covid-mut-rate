from Bio import SeqIO
from Bio import Entrez
Entrez.email = "jhwnkim@mit.edu"  # Always tell NCBI who you are

import traceback

# Read Covid19 reference sequence
ref = SeqIO.read("./data/ref_sequence.gb", "genbank")
print('Reference Covid sequence')
print(ref.id)
print(repr(ref.seq))
print(len(ref.seq))


def get_meta(record_id='NC_045512.2'): # default is the reference
    handle = Entrez.efetch(db="nuccore", id=record_id, rettype="gb", retmode="xml")
    fetch_record = Entrez.read(handle)
    handle.close()

    metadata = {
        'id': record_id,
        'submit-date': fetch_record[0]['GBSeq_create-date'],
    }

    for qualifier in fetch_record[0]['GBSeq_feature-table'][0]['GBFeature_quals']:
        if qualifier['GBQualifier_name'] == 'country':
            metadata['country'] = qualifier['GBQualifier_value']
        if qualifier['GBQualifier_name'] == 'collection_date':
            metadata['collect-date'] = qualifier['GBQualifier_value']

    print(metadata)
    return metadata


from Bio import Align
from Bio.Align.substitution_matrices import Array

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
    mut_rate = frequency / len(seq1)
    print(mut_rate)

    return mut_rate



# Read downloaded sequence file from NCBI GenBank Virus site
# records = list( SeqIO.parse("./data/MA-sequences.fasta", "fasta") )
# records = list( SeqIO.parse("./data/MA-sequences-1-toy.fasta", "fasta") )
records = list( SeqIO.parse("./data/MA-sequences-1.fasta", "fasta") )

metadata = []
mutarray = []

import time
start = time.time()
for idx, record in enumerate(records):
    print('\n{} of {} records'.format(idx+1, len(records)))
    try:
        meta = get_meta(record.id)
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

# Divide mutation rate by counts
for idx, idlist in enumerate(ids):
    mutarray_avg[idx] = mutarray_avg[idx]/len(idlist)

    print(dates[idx])
    print(idlist)
    print(mutarray_avg[idx])


# print('{:e}'.format(mutarray_avg[0]['C', 'T']))
# Save to file
import pickle

outfile = './data/MA-sequences-1.dat'
data = {
    'dates': dates,
    'idlist': idlist,
    'mutarray_avg': mutarray_avg
}
with open(outfile, "wb") as f:
    pickle.dump(data, f)


print('Run time took {}'.format(time.strftime("%H:%M:%S", time.gmtime(time.time()-start))))
