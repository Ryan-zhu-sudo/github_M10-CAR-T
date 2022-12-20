# amino acid ratios at each position of the statistical sequence

import os
import sys
from pathlib import Path
from Bio import SeqIO
from collections import Counter
import pandas as pd
import numpy as np
from itertools import chain


def is_low_freq_seq(seq, filter_aminos, max_lows=3):
    """
    seq: YKKY
    filter_aminos: {0: {Y, K}, 1: {}, }
    """
    fs = 0
    for n, a in enumerate(seq):
        filtas = filter_aminos[n]
        if a in filtas:
            fs += 1
    if fs > max_lows:
        return True
    return False


def get_seq_map(fp):
    res = {}
    fa = SeqIO.parse(fp, 'fasta')
    for recd in fa:
        name = recd.name
        seq = str(recd.seq)
        res[name] = seq
    return res


def mk_group(fp):
    grp_map = {}
    with open(fp) as fr:
        for line in fr:
            lst = line.strip().split('\t')
            grp_map[lst[0]] = lst[1]


def fetch_seqs(dirin, start_pos=0, end_pos=10):
    seqs = []
    aminos = set()
    files = chain(Path(dirin).rglob('*.txt'), Path(dirin).rglob('*.fas'))
    for f in files:
        print("process file: ", f)
        fa = SeqIO.parse(f, 'fasta')
        for recd in fa:
            # name = recd.name
            seq = str(recd.seq)
            selected_seq = seq[start_pos:end_pos+1]
            seqs.append(selected_seq)
            aminos.update(set(selected_seq))
    return seqs, list(aminos)


def fetch_files(files, start_pos=0, end_pos=10):
    seqs = []
    aminos = set()
    for f in files:
        print("process file: ", f)
        fa = SeqIO.parse(f, 'fasta')
        for recd in fa:
            # name = recd.name
            seq = str(recd.seq)
            selected_seq = seq[start_pos:end_pos+1]
            seqs.append(selected_seq)
            aminos.update(set(selected_seq))
    return seqs, list(aminos)


def fetch_oneseq(fpath, start_pos=0, end_pos=10):
    seqs = []
    aminos = set()
    print("process file: ", fpath)
    fa = SeqIO.parse(fpath, 'fasta')
    for recd in fa:
        seq = str(recd.seq)
        selected_seq = seq[start_pos:end_pos+1]
        seqs.append(selected_seq)
        aminos.update(set(selected_seq))
    return seqs, list(aminos)


def get_count(seqs, aminos, output=None,
              start_pos=0, end_pos=10,
              add_major=True):
    length = end_pos - start_pos + 1
    rows = []
    for i in range(length):
        ams = []
        for sq in seqs:
            if i <= len(sq)-1:
                ams.append(sq[i])
        _ct = dict(Counter(ams))
        row = []
        # _ct = dict(Counter([j[i] for j in seqs]))

        for a in aminos:
            row.append(_ct.get(a, 0))

        if add_major:
            _sum = sum(row)
            _max = max(row)
            a = aminos[row.index(_max)]
            # row.append(_sum)
            # row.append(_max/_sum)
            row.append(f'{a}|{round(_max/_sum,2)}')

        rows.append(row)

    if output is not None:
        df = pd.DataFrame(np.array(rows).T)
        df.index = aminos+['major']
        df.columns = range(start_pos, end_pos+1)
        df.to_csv(output)
        print('output: {}'.format(output))
    return rows


def run_stat(fpath, start_pos=0, end_pos=10, add_major=True):
    seqs, aminos = fetch_oneseq(fpath, start_pos=start_pos, end_pos=end_pos)
    rows = get_count(seqs, aminos, start_pos=start_pos, end_pos=end_pos, add_major=add_major)
    return rows, aminos


def run_stats(files, start_pos=0, end_pos=10, add_major=True):
    seqs, aminos = fetch_files(files, start_pos=start_pos, end_pos=end_pos)
    rows = get_count(seqs, aminos, start_pos=start_pos, end_pos=end_pos, add_major=add_major)
    return rows, aminos


if __name__ == "__main__":
    # 位置 0 起始
    fdir = sys.argv[1]
    START_POS = int(sys.argv[2])
    END_POS = int(sys.argv[3])
    outdir = sys.argv[4]

    # fdir = 'tests/PFB004_AD/before'
    # # fdir = r'H:\Projects\MyProjects\MYY_20220804\plotMotif\Input_file\1th\AD\before'
    # START_POS = 0
    # END_POS = 346
    # outdir = 'tests/PFB004_AD'

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    dr = os.path.basename(fdir.rstrip('/'))
    output_ = os.path.join(outdir, f'{dr}-seqstats.csv')
    seqs_, aminos_ = fetch_seqs(fdir, start_pos=START_POS, end_pos=END_POS)
    get_count(seqs_, aminos_, output_, start_pos=START_POS, end_pos=END_POS)
