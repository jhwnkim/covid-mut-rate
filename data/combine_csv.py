import pandas as pd
import sys
import glob

try:
    if len(sys.argv) > 1:
        directory = sys.argv[1]
        if directory[-1] != '/':
            directory = directory + '/'
    else:
        directory = './'
    csvfiles = glob.glob(directory+'*.csv')
except:
    print('csv results not available')
    exit()
else:
    print('Combining:')
    print(csvfiles)

print('\nReading {}'.format(csvfiles[0]))
allcsv = pd.read_csv(csvfiles[0])

for csvfile in csvfiles[1:]:
    print('Reading {}'.format(csvfile))
    allcsv = pd.concat([allcsv, pd.read_csv(csvfile)])

print('Sorting')
allcsv = allcsv.sort_values('Dates')

outfile = csvfiles[0][:-5]+'all.csv'
print('Writing to {}'.format(outfile))
allcsv.to_csv(outfile, index=False)
