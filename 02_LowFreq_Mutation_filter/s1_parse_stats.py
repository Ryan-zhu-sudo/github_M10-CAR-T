import pandas as pd
import os

from p1_stat_seq import fetch_seqs, get_count, fetch_oneseq, run_stat, run_stats


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


seq_lens = {'AD': 347, 'ID': 310}


def compare_two_seqs(x, y):
    diff = []
    for n, i in enumerate(x):
        allele, freq = i.split('|')
        y_c, y_freq = y[n].split('|')
        if allele != y_c:
            diff.append(f'{n}:{i}:{y[n]}')
    return diff


def get_ctrl_feq(diff, bf_freq, af_freqs):
    cts = []
    for m in diff:
        pos, before, after = m.split(':')
        # bf_c, bf_freq = before.split('|')
        # af_c, af_freq = after.split('|')
        ct_bf = bf_freq[int(pos)]
        ct_af = af_freqs[int(pos)]

        cb = ct_bf.split('|')[0]
        ca = ct_af.split('|')[0]
        tag = 'same'
        if ca != cb:
            tag = 'diff'

        cts.append(f'{pos}:{ct_bf}:{ct_af}:{tag}')
    return cts

# 对照组
control_data = {}
fp = 'control/results.xlsx'
dfs = pd.read_excel(fp, sheet_name=None)
sts = sorted(dfs.keys())
for sheet in sts:
    df = dfs[sheet]
    if sheet != 'Sheet':
        print(sheet)
        if 'AD' in sheet:
            tp = 'AD'
        else:
            tp = 'ID'
        if tp not in control_data:
            control_data[tp] = {}
        df['stage'] = df['sample'].apply(lambda x: status_map[x.split('_')[-1]])
        before_files = df['cluster_seq_file'][df['stage'] == 'Before'].tolist()
        after_files = df['cluster_seq_file'][df['stage'] == 'After'].tolist()
        before_files_ = [i.replace('3th-new/3th-cluster0.95-stats', 'control') for i in before_files]
        after_files_ = [i.replace('3th-new/3th-cluster0.95-stats', 'control') for i in after_files]

        before_srows = run_stats(before_files_, start_pos=0, end_pos=seq_lens[tp]-1)
        after_srows = run_stats(after_files_, start_pos=0, end_pos=seq_lens[tp]-1)
        before_major_frqs = [i[-1] for i in before_srows]
        after_major_frqs = [i[-1] for i in after_srows]
        control_data[tp][sheet] = [before_major_frqs, after_major_frqs]
print(control_data)

######################################################################
dns = ['1th-new/1th-cluster0.95-stats', '3th-new/3th-cluster0.95-stats']

output = './mutations.csv'
fw = open(output, 'w', encoding='utf-8')
header = [
    'batch', 'sample', 'group', 'format', 'mutations', 'control1', 'control2', 'control3',
]
fw.write(','.join(header) + '\n')

for dn in dns:
    fp = dn + '/results.xlsx'

    dfs = pd.read_excel(fp, sheet_name=None)
    for sheet, df in dfs.items():
        if sheet != 'Sheet':
            print(sheet)
            if 'AD' in sheet:
                tp = 'AD'
            else:
                tp = 'ID'
            # print(df.head())
            before_max = df['before_count'].max()
            after_max = df['after_count'].max()
            dfa = df[(df['before_count'] == before_max) | (df['after_count'] == after_max)].copy()
            dfa['stage'] = df['sample'].apply(lambda x: status_map[x.split('_')[-1]])
            before_seq = dfa['seq'][dfa['stage'] == 'Before'].squeeze()
            after_seq = dfa['seq'][dfa['stage'] == 'After'].squeeze()
            before_file = dfa['cluster_seq_file'][dfa['stage'] == 'Before'].squeeze()
            after_file = dfa['cluster_seq_file'][dfa['stage'] == 'After'].squeeze()

            before_srows = run_stat(before_file, start_pos=0, end_pos=seq_lens[tp]-1)
            after_srows = run_stat(after_file, start_pos=0, end_pos=seq_lens[tp]-1)
            before_major_frqs = [i[-1] for i in before_srows]
            after_major_frqs = [i[-1] for i in after_srows]
            # print(before_major_frqs)
            # print(after_major_frqs)

            diffs = compare_two_seqs(before_major_frqs, after_major_frqs)
            # print(diffs)
            control_d = control_data[tp]

            row = [dn, sheet, tp, 'pos(0base):Before:After', ";".join(diffs)]
            for k, v in control_d.items():
                kct = get_ctrl_feq(diffs, v[0], v[1])
                row.append(k + "@" + ";".join(kct))

            # fw.write(f'{dn}\t{sheet}\t{tp}\tpos(0base):Before:After\t{";".join(diffs)}\n')
            fw.write(','.join(row) + '\n')
fw.close()



