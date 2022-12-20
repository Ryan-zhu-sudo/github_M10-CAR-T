import os.path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from collections import Counter

from b1_create_query import INTERVAL

# Allow blast difference length as a percentage of query length
BLAST_DIFF_RATE = 0.5

# blast minimum percentage
MIN_BLAST_PERCENT = 90.

# PCR regions, reference HIV-genome.docx
AD_START = 2244
AD_END = 3314

ID_START = 4164
ID_END = 5093

LENS = {
    'AD': AD_END-AD_START+1,
    'ID': ID_END-ID_START+1,
}

# Parameters for the proportion of sequences intercepted at the beginning and end of the sequence, due to incomplete sequencing, some sequences do not start from the reference starting point, the list indicates the maximum proportion of incomplete sequenced sequences retained after interception filtering
# The first element indicates the starting position, the second indicates the end position, the higher the ratio means the longer the length of the cropped sequence, the shorter the final retained sequence
# The intercept length will eventually be corrected upwards to a multiple of 3
# This parameter is used for pre-testing and screening for suitable lengths
TRIMS = {
    'AD': [0.8, 0.5],  # start, end
    'ID': [0.8, 0.5]
}


TRIMS_LENS = {
    'AD': [0, 0],  # start, end
    'ID': [0, 0]
}


def get_fa_map(fp):
    res = {}
    fa = SeqIO.parse(fp, 'fasta')
    for recd in fa:
        name = recd.name
        seq = recd.seq
        res[name] = seq
    return res


def get_rev_compl_seq(seq):
    s = Seq(seq, IUPAC.unambiguous_dna)
    srvc = s.reverse_complement()
    return srvc._data


def cal_cut_rate(ct, rate=0.8):
    _n = sum(ct.values())
    lim = _n*rate
    sm = 0
    vs = []
    for i in sorted(ct.keys()):
        if not sm > lim:
            sm += ct[i]
            vs.append(i)
        else:
            break
    cv = vs[-1]
    _m = cv % 3
    if _m != 0:
        cv = cv + (3-_m)
    return cv


if __name__ == "__main__":
    fa_dir = '../03_SeqCount'
    outdir = '../output'
    result_dir = '../result'
    if not os.path.exists(result_dir):
        os.mkdir(result_dir)

    files = [i for i in os.listdir(outdir) if i.endswith('.blast')]
    for f in files:
        print(f)

        datas = []
        start_diffs = []
        end_diffs = []

        sample = f.split('.')[0]
        tp = f.split('.')[1]

        target_length = LENS[tp]

        output = os.path.join(outdir, f'{sample}.{tp}.blast.parse.txt')
        fw = open(output, 'w')
        header = 'seq_id\tfa_id\tstrand\tstart_diff\tend_diff\tfa_len\tfa_seq\n'
        fw.write(header)

        # 0-base
        query_start, query_end = INTERVAL[tp]
        query_length = query_end - query_start + 1
        blast_length_diff = int(query_length * BLAST_DIFF_RATE)

        fp_blast = os.path.join(outdir, f)
        fp_fa = os.path.join(fa_dir, sample, f'{sample}.target.fa')
        fa_map = get_fa_map(fp_fa)

        n = 1
        with open(fp_blast) as fr:
            for line in fr:
                line = line.strip()
                lst = line.split('\t')
                fa_id = lst[1]
                blast_length = int(lst[3])
                blast_percent = float(lst[2])

                blast_query_start, blast_query_end = int(lst[6]), int(lst[7])
                blast_db_start, blast_db_end = int(lst[8]), int(lst[9])
                if abs(query_length - blast_length) < blast_length_diff and blast_percent >= MIN_BLAST_PERCENT:

                    offset = blast_query_start - 1

                    raw_seq = fa_map[fa_id]
                    raw_seq_length = len(raw_seq)

                    if blast_db_start < blast_db_end:
                        strand = "+"
                    else:
                        strand = "-"
                        raw_seq = get_rev_compl_seq(str(raw_seq))
                        blast_db_start = raw_seq_length - blast_db_start + 1
                        blast_db_end = raw_seq_length - blast_db_end + 1

                    _start = blast_db_start - offset

                    trim_start = _start
                    polish = ''
                    if _start <= 0:
                        trim_start = 1
                        for i in range(_start, 1):
                            polish += '-'
                    start_trim_seq = polish + raw_seq[trim_start-1:]

                    # print(n, fa_id, strand, start_trim_seq[:50])

                    trim_seq = start_trim_seq[:target_length]
                    seq_length = len(trim_seq)

                    if seq_length > target_length / 2:

                        start_diff = len(polish)
                        end_diff = target_length - seq_length

                        if start_diff > 0:
                            start_diffs.append(start_diff)
                        if end_diff > 0:
                            end_diffs.append(end_diff)

                        seq_id = '{}_{}'.format(sample, n)
                        row = [str(s) for s in [seq_id, fa_id, strand, start_diff, end_diff, seq_length, trim_seq]]

                        fw.write('\t'.join(row)+'\n')

                        datas.append(row)

                        n += 1
        fw.close()

        if datas:
            # Interception filtering, according to the overall proportion to determine the interception position, certain sequence start and end sequence is not measured, interception still exist after the gap of the sequence will be filtered
            output = os.path.join(result_dir, f'{sample}.{tp}-dna.fa')
            fw = open(output, 'w')

            # for test
            # start_cut_rate, end_cut_rate = TRIMS[tp]
            # stc = Counter(start_diffs)
            # edc = Counter(end_diffs)
            # if stc:
            #     start_cut_length = cal_cut_rate(stc, start_cut_rate)
            #     print('start_cut_length: ', start_cut_length)
            # if edc:
            #     end_cut_length = cal_cut_rate(edc, end_cut_rate)
            #     print('end_cut_length: ', end_cut_length)

            start_cut_length = TRIMS_LENS[tp][0]
            end_cut_length = TRIMS_LENS[tp][1]

            for rw in datas:
                if int(rw[3]) <= start_cut_length and int(rw[4]) <= end_cut_length:
                    seq_id, fa_id, strand = rw[:3]
                    seq_len = int(rw[-2])
                    if end_cut_length == int(rw[4]):
                        cut_seq = rw[-1][start_cut_length+1:]
                    else:
                        cut_seq = rw[-1][start_cut_length+1: 1-(end_cut_length-int(rw[4])+1)]
                    # print(rw[0], len(cut_seq))
                    new_row = []
                    fw.write(f'>{seq_id}\n')
                    fw.write(f'{cut_seq}\n')

            fw.close()
        else:
            print('no data!')


