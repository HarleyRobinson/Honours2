import csv
import re
f= open('lipidraftSig.csv', 'r')
new= []

for line in f:
    result= line.split(",")
    m= re.findall("GN=(\w+)", str(result))
    if len(m)==1:
        new.append(m)
with open('RaftSigParsed.csv', 'w') as f:
    writer=csv.writer(f, delimiter= ',')
    writer.writerows(new)
