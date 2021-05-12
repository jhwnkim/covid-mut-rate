import pandas as pd
import sys

print(sys.argv)
if len(sys.argv) > 1:
    infile = sys.argv[1]
else:
    infile = "./old/MA-sequences-2-toy.fasta"

print('Reading csv file {}'.format(infile))
df = pd.read_csv(infile)

print('Performing day averages')
row = 0
while row < len(df):
	date = df.iloc[row,0]
	location = df.iloc[row,2]

	dfday = df.loc[df['Dates']==date]
	# print(dfday)

	dfout_row = pd.DataFrame({'Dates': [date], 'Geolocation': [location]})

	dfout_row = pd.concat([dfout_row, dfday.iloc[:,3:].sum().to_frame().transpose()], axis=1)

	dfout_row['N'] = len(dfday)
	if row == 0:
		dfout = dfout_row
	else:
		dfout = pd.concat([dfout, dfout_row])
	row = row + len(dfday)

outfile = infile[:-4]+'-avg.csv'
print('Writing csv file {}'.format(outfile))

dfout.to_csv(outfile, index=False)
