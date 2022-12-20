# _*_ coding: utf-8 _*_
# @Time     : 2022/12/9 0009 14:04
# @Author   : Yang Liuke

import sys
from Bio import SeqIO

# from pathlib import Path
import os

fdir = sys.argv[1]
# fdir = '../tests/translate'

fdiro = fdir+'-filter'
if not os.path.exists(fdiro):
    os.mkdir(fdiro)

for f in os.listdir(fdir):
    fp = os.path.join(fdir, f)
    output = os.path.join(fdiro, f)
    fw = open(output, 'w')
    exist_ids = set()
    with open(fp) as fr:
        fa = SeqIO.parse(fr, 'fasta')
        for recd in fa:
            name = recd.name
            seq = recd.seq
            seq_id = name.rsplit('_', 1)[0]
            if '*' not in seq and seq_id not in exist_ids:
                fw.write('>{}\n{}\n'.format(name, seq))
                exist_ids.add(seq_id)

    fw.close()


