import csv
import re
f= open('ProteosomeInsig.csv', 'r')
new= []

for line in f:
    result= line.split(",")
    m= re.findall("GN=(\w+)", str(result))
    if len(m)==1:
        new.append(m)
with open('ExoinsigParsed.csv', 'w') as f:
    writer=csv.writer(f, delimiter= ',')
    writer.writerows(new)
