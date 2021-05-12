import pandas as pd

import glob

try:
    csvfiles = glob.glob('./*.csv')
except:
    print('csv results not available')
    exit()
else:
    print(csvfiles)

allcsv = pd.read_csv(csvfiles[0])

for file in csvfiles[1:]:
    allcsv = pd.concat([allcsv, pd.read_csv(file)])


allcsv = allcsv.sort_values('Dates')

allcsv.to_csv(csvfiles[0][:-5]+'all.csv', index=False)
