import csv
import re
import os
os.chdir('C:\\Users\\hrobinson\\Documents\\FromVirtualBox')
f= open('AllMatchedToMotifs.csv', 'r')
new= []
for line in f:
    result= line.replace(' ', ',')
    new.append(result)
new=(", ".join(new))
with open('MatchedMotifsAllMirs.csv', 'w') as f:
    writer=f.write(str(new))
    #writer.writerows(result)
