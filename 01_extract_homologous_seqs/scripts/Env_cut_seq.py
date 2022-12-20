import os.path
import re
from hashlib import md5
from pathlib import Path
from Bio import SeqIO


N = 25

LENGTH = 2568
# TARGET = "ATGAGAGTGA"
TARGET = "ATGAGA"


def get_md5(x):
    mv = md5()
    mv.update(x.encode())
    return mv.hexdigest()


def get_sample_fa(path, outdir):
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    files = Path(path).rglob('*target.fa')
    for f in files:
        # print("process: {}".format(f))
        tp = f.name.split('.')[0]
        output = os.path.join(outdir, f'target-{tp}.fa')

        fwn = open(output, 'w')
        with open(output, 'w') as fw:
            fa = SeqIO.parse(f, 'fasta')
            n = 0
            m = 0
            d = 0
            eseqs = set()
            for recd in fa:
                name = recd.name
                seq = str(recd.seq)
                if len(seq) > 3000:
                    try:
                        # h, t = re.split(re.compile(TARGET), seq)
                        h, *ss = seq.split(TARGET)
                        n += 1
                        nseq = (TARGET+''.join(ss))[:LENGTH]
                        # print(len(nseq))
                        # print(nseq)
                        hlen = len(h)
                        # print(hlen)
                        # if not 250 < hlen < 310:
                        #     print(f, name, nseq[:50], sep='\t')
                        m5 = get_md5(nseq)
                        if len(nseq) == LENGTH and 270 < hlen < 300:
                            if m5 not in eseqs:
                                nname = '{}_{}'.format(tp, n)
                                fwn.write(f'>{nname}\n{nseq}\n')
                                # print(f, name, nseq[:10], sep='\t')
#                                eseqs.add(m5)
                            else:
                                print(f, name, nseq[:10], sep='\t')
                                d += 1
                    except:
                        # print(f, name)
                        m += 1
                    # fwn.write(f'>{name}\n{seq}\n')
        fwn.close()
        print(n, m, d)
        # print("output: {}".format(output))


if __name__ == "__main__":
    fdir_ = '../03_Env_SeqCount'
    outdir_ = '../Env_output'
    get_sample_fa(fdir_, outdir_)

