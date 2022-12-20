INTERVAL = {
    'AD': [1, 50],
    'ID': [1, 50],
}


for tp in ['AD', 'ID']:
    fa = '../query/{}.fa'.format(tp)
    with open(fa) as fr:
        _id = fr.readline().split(' ')[0][1:]
        _sq = ''.join(fr.readlines()).replace('\n', '').strip()
        start_, end = INTERVAL[tp]
        start = start_ - 1
        output = '../query/{}-query.fa'.format(tp)
        with open(output, 'w') as fw:
            fw.write('>{}\n'.format(tp))
            fw.write(_sq[start:end]+'\n')