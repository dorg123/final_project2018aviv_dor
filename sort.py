import fasta
reader = fasta.TabReader('pairs_resolution.txt')
data = sorted(reader.data, key=lambda x: float(x['res']), reverse=True)
with open('pairs_by_resolution.txt', 'w') as f:
    f.write('subseq1\tsubseq2\tseqnum\thist400\tres\tavlen\n')
    for d in data:
        f.write('{d[subseq1]}\t{d[subseq2]}\t{d[seqnum]}\t{d[hist400]}\t{d[res]}\t{d[avlen]}\n'.format(d=d))
print('done')
