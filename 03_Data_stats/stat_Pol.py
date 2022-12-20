import os
from pathlib import Path
from Bio import SeqIO
from hashlib import md5
import sys
from openpyxl import Workbook
from collections import Counter

# dir_clstr = '1th-cluster'
# dir_fasta = '1th-filter'

dir_clstr = sys.argv[1]
dir_fasta = sys.argv[2]
#

outdir = dir_clstr+'-stats'
if not os.path.exists(outdir):
    os.mkdir(outdir)

seq_dir = os.path.join(outdir, 'seqs')
if not os.path.exists(seq_dir):
    os.mkdir(seq_dir)


print("##0")
# cd hit
#################################################
files = [i for i in os.listdir(dir_clstr) if i.endswith('.clstr')]
datas = []
cluster_counts = {}
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
            # print(_i, _p)
            ck = f"{f}#{key}"
            if ck not in cluster_counts:
                cluster_counts[ck] = 0
            cluster_counts[ck] += 1
            datas.append((_i, [f, _p, key]))

cdhit_map = {}
for k, v in datas:
    cdhit_map[k] = v

print("##1")
#################################################
output = os.path.join(outdir, 'pro_seqs.csv')
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
fw = open(output, 'w')
fw.write('sample,group,status,type,seq_id,seq_md5,cdhit_file,cdhit_percent,cdhit_cluster,cluster_count,seq\n')

files = Path(dir_fasta).rglob('*.txt')
sdatas = {}
for f in files:
    print(f)
    fname = f.name
    fs = fname.split('.')[0]
    sample, typ = fs.rsplit('_', 1)

    sg = sample.split('_')[0]
    tg = sample.split('_')[-1]

    group = f"{sg}_{typ}"
    if sample in ['AD3', 'AD24']:
        group = 'AD3_AD24'
    elif sample in ['ID4', 'ID25']:
        group = 'ID4_ID25'

    if group not in sdatas:
        sdatas[group] = []

    # print(group, typ, tg)

    fa = SeqIO.parse(f, 'fasta')
    n_seqs = 0
    for recd in fa:
        n_seqs += 1
        m = md5()
        m.update(str(recd.seq).encode())
        mv = m.hexdigest()
        seq_id = recd.id
        seq = str(recd.seq)
        status = status_map[tg]
        cdhit_file,cdhit_percent,cdhit_cluster = cdhit_map[seq_id]

        ck = f"{cdhit_file}#{cdhit_cluster}"
        cluster_count = cluster_counts[ck]
        row = [
            sample, group, status, typ,
            seq_id, mv,
            cdhit_file, cdhit_percent, cdhit_cluster,
            cluster_count, seq
        ]
        sdatas[group].append(row)
        # print(*row, sep=',', file=fw)
fw.close()


print("##2")
###############################################
# cdhit_cluster	before_count	after_count	sample	seq
all_seqs = {}
cluster_seqs = {}
cluster_seqs_seq = {}
res = {}
for rg, rows in sdatas.items():
    print(rg)
    # print(rows)

    if rg not in res:
        res[rg] = {}

    if rg not in cluster_seqs:
        cluster_seqs[rg] = {}

    if rg not in all_seqs:
        all_seqs[rg] = {'Before': [], 'After': []}

    for row in rows:
        ck = f"{row[-5]}#{row[-3]}"
        print(ck)
        if ck not in res[rg]:
            cdhit_cluster = ck.split('#')[1]
            res[rg][ck] = {
                    'cdhit_cluster': cdhit_cluster,
                    'before_count': 0,
                    'after_count': 0,
                    'sample': '',
                    'seq': '',
                    'status': ''
                    }
            cluster_seqs[ck] = []
            cluster_seqs_seq[ck] = []

        (sample, group, status, typ,
         seq_id, mv,
         cdhit_file, cdhit_percent, cdhit_cluster,
         cluster_count, seq) = row

        cluster_seqs[ck].append(mv)
        cluster_seqs_seq[ck].append((seq_id, seq))

        if status == 'Before':
            res[rg][ck]['before_count'] += 1
            all_seqs[rg]['Before'].append(mv)

        elif status == 'After':
            res[rg][ck]['after_count'] += 1
            all_seqs[rg]['After'].append(mv)

        if cdhit_percent == '*':
            res[rg][ck]['seq'] = seq

        if not res[rg][ck]['sample']:
            res[rg][ck]['sample'] = sample

        if not res[rg][ck]['status']:
            res[rg][ck]['status'] = status

# print(res)
# print(cluster_seqs_seq.keys())

#################################################


def get_intersect(x1, x2):
    _ct = Counter(x2)
    common = set(x1) & set(x2)
    if common:
        s = sum([_ct[i] for i in common])
        # print(common, s)
        return s
    else:
        return 0


print("##3")
# #################################################
output = os.path.join(outdir, 'results.xlsx')
wb = Workbook()
_n = 0
for k, v in res.items():
    # print(k)

    seq_outdir = os.path.join(seq_dir, k)
    if not os.path.exists(seq_outdir):
        os.mkdir(seq_outdir)

    represents_fa = os.path.join(seq_dir, 'represents.{}.fas'.format(k))
    fwr = open(represents_fa, 'w')

    st = wb.create_sheet(k, _n)
    _n += 1
    header = 'cdhit_cluster,before_count,after_count,sample,seq,cluster_seq_file'.split(',')
    st.append(header)
    before_seqs = all_seqs[k]['Before']
    after_seqs = all_seqs[k]['After']
    # print(len(before_seqs), len(after_seqs))
    for _k, _v in v.items():
        print(_k)

        cluster_fa = os.path.join(seq_outdir, f'{_k}.fas')
        seq_tup = cluster_seqs_seq[_k]
        # print(len(seq_tup))
        with open(cluster_fa, 'w') as fwcc:
            for sk, sv in seq_tup:
                fwcc.write(f'>{sk}\n')
                fwcc.write(f'{sv}\n')

        _vdata = list(_v.values())
        status = _vdata[-1]
        _k_seqs = cluster_seqs[_k]
        # print(status, len(_k_seqs))

        rw = _vdata[:-1]

        if status == 'Before':
            after_count = get_intersect(_k_seqs, after_seqs)
            rw[2] = after_count

        elif status == 'After':
            before_count = get_intersect(_k_seqs, before_seqs)
            rw[1] = before_count

        rw.append(cluster_fa.replace('../', ''))
        st.append(rw)

        rseq_id = '{}_{}'.format(rw[3], rw[0].replace(' ', ''))
        rseq = rw[-2]
        fwr.write('>{}\n'.format(rseq_id))
        fwr.write('{}\n'.format(rseq))

    fwr.close()
wb.save(output)

if os.path.exists(output):
    print('output: {}'.format(output))

# #################################################
# df = pd.read_csv(output)
# for i in ['AD', 'ID']:
#     output2 = f'../{i}-intersect.csv'
#     fw2 = open(output2, 'w')
#     dfi = df[df['group'] == i]
#     grp = dfi.groupby(['seq_md5'])
#     dfs = []
#     for g, dfg in grp:
#         df_b = dfg[dfg['status'] == 'Before']
#         df_a = dfg[dfg['status'] == 'After']
#         if df_b.shape[0] > 0 and df_a.shape[0] > 0:
#             # print(g, dfg.shape)
#             dfs.append(dfg)
#     dfm = pd.concat(dfs)
#     dfm.to_csv(output2, index=False)
# fw2.close()

