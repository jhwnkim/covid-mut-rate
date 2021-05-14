import pandas as pd
import sys

print(sys.argv)
if len(sys.argv) > 1:
    infiles = [sys.argv[1]]
else:
    infiles = ["AL-sequences-all.csv", \
				"AZ-sequences-all.csv", \
				"CA-sequences-all.csv", \
				"CT-sequences-all.csv", \
				"FL-sequences-all.csv", \
				"GA-sequences-all.csv", \
				"IL-sequences-all.csv", \
				"MA-sequences-all.csv", \
				"MD-sequences-all.csv", \
				"MN-sequences-all.csv", \
				"NC-sequences-all.csv", \
				"NJ-sequences-all.csv", \
				"NM-sequences-all.csv", \
				"NY-sequences-all.csv", \
				"OH-sequences-all.csv", \
				"PA-sequences-all.csv", \
				"RI-sequences-all.csv", \
				"SC-sequences-all.csv", \
				"TN-sequences-all.csv", \
				"TX-sequences-all.csv", \
				"VA-sequences-all.csv", \
				"WA-sequences-all.csv", \
				"WI-sequences-all.csv", ]

for infile in infiles:
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

		dfout_row = pd.concat([dfout_row, dfday.iloc[:,3:].sum().to_frame().transpose()/len(dfday)], axis=1)

		dfout_row['N'] = len(dfday)
		if row == 0:
			dfout = dfout_row
		else:
			dfout = pd.concat([dfout, dfout_row])
		row = row + len(dfday)

	outfile = infile[:-4]+'-avg.csv'
	print('Writing csv file {}'.format(outfile))

	dfout.to_csv(outfile, index=False)
