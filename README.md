# covid-mut-rate
Covid19 Mutation Rate Prediction Project

## Environment Setup

```
$ conda create --name cmr
$ conda activate cmr
$ conda install jupyter ipython pandas
$ conda install -c conda-forge biopython
```

## Datasets

* ./data/MA-sequences-1.fasta: Fasta list of following search

1. Virus: Sever acute respiratory syndrome coronavirus 2(SARS-Cov-2), taxid:2697049

2. Sequence length: 29161-29903

3. Nucleotide completeness: completeness

4. Geographic region: USA: MA

5. Collection Date: Apr 26,2021 - May 10,2021

Fasta definition line: Accession, Genbank Title


* ./data/MA-sequences-2.fasta: Fasta list of following search

1. Virus: Sever acute respiratory syndrome coronavirus 2(SARS-Cov-2), taxid:2697049

2. Sequence length: 29161-29903

3. Nucleotide completeness: completeness

4. Geographic region: USA: MA

5. Collection Date: Apr 10,2021 - May 10,2021

Fasta definition line: Accession, Genbank Title, Geo Location, Collection Date
