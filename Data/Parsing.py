#import os
#os.chdir(directory)
import csv
f= open('HumanProtInter.csv', 'r')
FUS= []
for line in f:
    result= line.split(",")
    if "HNRNPK" == result[7]:
        FUS.append(result)
    if "HNRNPK"== result[8]:
        FUS.append(result)

with open('HNRNPKinter.csv', 'w') as f:
    writer=csv.writer(f, delimiter=',')
    writer.writerows(FUS)

