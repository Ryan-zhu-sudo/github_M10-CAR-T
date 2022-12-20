# _*_ coding: utf-8 _*_
# @Time     : 2022/9/20 0020 9:08
# @Author   : yangliuke

import os
import sys
# import re
from Bio import SeqIO
from openpyxl import Workbook

# dirs = [('1th-combine-cluster', '1th-filter'),
#         ('3th-combine-cluster', '3th-filter'),
#         ('env-combine-cluster', 'env-filter')
#         ]

dir_clstr = sys.argv[1]
dir_fasta = sys.argv[2]

# dir_clstr = '1th-combine-cluster'
# dir_fasta = '1th-filter'

mdir = dir_clstr + '-stats'
if not os.path.exists(mdir):
    os.mkdir(mdir)

seq_dir = os.path.join(mdir, 'seqs')
if not os.path.exists(seq_dir):
    os.mkdir(seq_dir)

output = os.path.join(mdir, 'result.xlsx')


#########################################
AFTER = 'After'
BEFORE = 'Before'

status_map = {
    'AD3': 'After',
    'AD24': 'Before',
    'D-2': 'Before',
    'D-3': 'Before',
    '2D0': 'Before',
    'shai': 'Before',
    'SHAI': 'Before',
    'after': 'After',
    'AFTER': 'After',
    'M4': 'After',
    'M5': 'After',
    'M2': 'After',
    'ID4': 'After',
    '21-7': 'After',
    '21-10': 'After',
    'ID25': 'Before',
    'D1': 'Before',
    '21-3': 'Before',
    '21-4': 'Before',
    'D-30': 'Before',
}


def get_sample_tp(x):
    ss = x.split('_')[0].split('-')
    if len(ss) == 1:
        sp = tp = ss[0]
    else:
        sp = ss[0]
        tp = '-'.join(ss[1:])
    return sp, tp


def _count(x):
    bc = 0
    ac = 0
    for xi in x:
        sp, tp = get_sample_tp(xi)
        st = status_map[tp]
        if st == AFTER:
            ac += 1
        else:
            bc += 1
    return bc, ac


def get_all_seqs(dpath):
    seqs = {}
    files_ = [i for i in os.listdir(dpath) if i.endswith('.txt')]
    for f_ in files_:
        fp_ = os.path.join(dpath, f_)
        with open(fp_) as fr_:
            fa_ = SeqIO.parse(fr_, 'fasta')
            for recd in fa_:
                seqs[recd.name] = str(recd.seq)
    return seqs


files = [i for i in os.listdir(dir_clstr) if i.endswith('.clstr')]
cluster_seqs = {}
key = None
for f in files:
    print(f)
    fr = open(os.path.join(dir_clstr, f))
    for line in fr:
        if line.startswith('>'):
            key = line.strip('\n>')
        else:
            n, row = line.strip().split('\t')
            _i = row.split('.')[0].split('>')[1]
            _p = row.split(' ')[-1]
            ck = f"{f}#{key}"

            if f not in cluster_seqs:
                cluster_seqs[f] = {}
            if key not in cluster_seqs[f]:
                cluster_seqs[f][key] = []
            cluster_seqs[f][key].append(_i)

#
all_seqs = get_all_seqs(dir_fasta)
wb = Workbook()
n = 0
for fn, cs in cluster_seqs.items():
    st_name = fn.split('_com')[0]
    st = wb.create_sheet(st_name, n)
    header = 'cdhit_cluster,before_count,after_count,sample,seq,cluster_seq_file'.split(',')
    st.append(header)

    fn_seq = os.path.join(seq_dir, 'represents.{}.fas'.format(st_name))
    fw = open(fn_seq, 'w')

    n += 1
    for c, cl in cs.items():
        before_count, after_count = _count(cl)
        represent_seq = all_seqs[cl[0]]
        fw.write(f'>{c}\n{represent_seq}\n')
        sample = st_name
        row = [c, before_count, after_count, sample, represent_seq, fn]
        st.append(row)

        fc_dir = os.path.join(seq_dir, st_name)
        if not os.path.exists(fc_dir):
            os.mkdir(fc_dir)
        fcout = os.path.join(fc_dir, '{}#{}.fas'.format(fn, c))
        with open(fcout, 'w') as fww:
            for sid in cl:
                fww.write(">{}\n{}\n".format(sid, all_seqs[sid]))

    fw.close()
wb.save(output)

if os.path.exists(output):
    print('output: {}'.format(output))

