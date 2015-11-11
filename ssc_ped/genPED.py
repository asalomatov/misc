#import pandas as pd
#xl_file = pd.ExcelFile('Eichler_table_S1.xlsx')
#mdf = xl_file.parse(xl_file.sheet_names[0])

import sys

sep = sys.argv[-1]
prob = []
sib = []
for b in [1,2,3,4,5]:
    prob.append('p'+str(b))
    sib.append('s'+str(b))
with open('SSC.ped', 'w') as fout:
    with open('Eichler_table_S1.csv', 'r') as fin:
        next(fin)
        for l in fin:
            sl = l.split(',')
            faID = sl[0]+sep+'fa'
            moID = sl[0]+sep+'mo'
            fa = '\t'.join([sl[0], faID, str(0), str(0), str(1), str(1), 'SSC'])
            mo = '\t'.join([sl[0], moID, str(0), str(0), str(2), str(1), 'SSC'])
            fout.write(fa+'\n')
            fout.write(mo+'\n')
            for p in prob:
                if p in l: 
                    sex = int(sl[2]=='M')
                    if sex == 0: sex = 2
                    proband = '\t'.join([sl[0], sl[0]+sep+p, faID, moID, str(sex), str(2), 'SSC'])
                    fout.write(proband+'\n')
            for s in sib:
                if s in l: 
                    sex = int(sl[2]=='M')
                    if sex == 0: sex = 2
                    sibling = '\t'.join([sl[0], sl[0]+sep+s, faID, moID, str(sex), str(1), 'SSC'])
                    fout.write(sibling+'\n')
