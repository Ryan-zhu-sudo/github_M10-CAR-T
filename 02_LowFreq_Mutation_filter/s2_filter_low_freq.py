# Filter sequences containing low frequency loci
import numpy as np
import pandas as pd
import os
import sys

from p1_stat_seq import (get_count, fetch_oneseq,
                         get_seq_map, is_low_freq_seq)

# fdir = './test_before'
# fdir = './1th-new/1th-cluster0.95-stats/seqs\AD3_AD24'
# fdir = './AD_rename-filter'
# start = 0
# end = 344
# min_freq = 0.01
# max_lows = 3

fdir = sys.argv[1]
start = int(sys.argv[2])
end = int(sys.argv[3])
min_freq = float(sys.argv[4])
max_lows = int(sys.argv[5])

outdir = fdir.rstrip('/').rstrip('\\') + '_filter_low_freq'
if not os.path.exists(outdir):
    os.mkdir(outdir)

files = [os.path.join(fdir, i) for i in os.listdir(fdir)]
for f in files:
    seqs, aminos = fetch_oneseq(f, start_pos=start, end_pos=end)
    rows = get_count(seqs, aminos, start_pos=start, end_pos=end, add_major=False)

    df = pd.DataFrame(np.array(rows).T)
    df.index = aminos
    df.to_excel(f'./output/{os.path.basename(f)}.xlsx')

    # print(aminos)
    # print(rows)

    filter_aminos = {}
    for n, row in enumerate(rows):
        if n not in filter_aminos:
            filter_aminos[n] = set()
        for m, v in enumerate(row):
            m_freq = v / sum(row)
            if 0 < m_freq < min_freq:
                filter_aminos[n].add(aminos[m])
    # print(filter_aminos)

    fn = os.path.basename(f).rsplit('.', 1)[0]
    pass_dir = os.path.join(outdir, 'passed')
    filter_dir = os.path.join(outdir, 'filtered')
    if not os.path.exists(pass_dir):
        os.mkdir(pass_dir)
    if not os.path.exists(filter_dir):
        os.mkdir(filter_dir)
    output_passed = os.path.join(pass_dir, f'{fn}.fas')
    output_filter = os.path.join(filter_dir, f'{fn}.fas')

    fw_pass = open(output_passed, 'w')
    fw_filter = open(output_filter, 'w')
    seq_map = get_seq_map(f)
    filtered_map = {}
    for nm, seq in seq_map.items():
        if not is_low_freq_seq(seq, filter_aminos, max_lows=max_lows):
            filtered_map[nm] = seq
            fw_pass.write(f">{nm}\n{seq}\n")
        else:
            fw_filter.write(f">{nm}\n{seq}\n")
    print(f"raw seqs：{len(seq_map)}, seqs after filtering：{len(filtered_map)}")
