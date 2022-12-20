import pandas as pd
import sys
import os

fin1 = sys.argv[1]
fin2 = sys.argv[2]

# fin1 = 'output/test_before-seqstats.csv'
# fin2 = 'output/test_after-seqstats.csv'

outdir = os.path.dirname(fin1)
output = os.path.join(outdir, 'amino_compare.csv')

df1 = pd.read_csv(fin1, index_col=0)
df2 = pd.read_csv(fin2, index_col=0)

d1 = dict(zip(df1.columns, df1.loc['major']))
d2 = dict(zip(df2.columns, df2.loc['major']))

with open(output, 'w', encoding='utf-8') as fw:
    fw.write(f'pos,p1_major,p2_major,p1_freq,p2_freq,type,diff_freq\n')
    for k, v in d1.items():
        v2 = d2.get(k)
        a1, f1 = v.split('|')
        a2, f2 = v2.split('|')
        diff_freq = float(f2) - float(f1)
        tp = 'same'
        if a1 != a2:
            tp = 'different'
        fw.write(f'{k},{a1},{a2},{f1},{f2},{tp},{diff_freq}\n')

print('output: {}'.format(output))
