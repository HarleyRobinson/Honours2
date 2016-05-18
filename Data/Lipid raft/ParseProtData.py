import csv
import re
import os
os.chdir('C:\\Users\\hrobinson\\Documents\\FromVirtualBox')
f= open('AllMirsInDataSet.fa', 'r')
new= []
for line in f:
    result= line.split("\t")
    new.append(result)
    
with open('ParsedMirs.fa', 'w') as f:
    writer=f.write(str(new))
    #writer.writerows(result)
