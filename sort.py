import fasta
# read the pairs file
reader = fasta.TabReader('pairs_resolution.txt')
data = sorted(reader.data, key=lambda x: float(x['res']), reverse=True)  # sort the data by resolution
with open('pairs_by_resolution20.txt', 'w') as f:
    f.write('subseq1\tsubseq2\tseqnum\thist400\tres\tavlen\n')
    for d in data:
        if len(d['subseq1']) == 20 and len(d['subseq2']) == 20:  # filter by primer length 20
            # and write it down with the real resolution
            f.write('{d[subseq1]}\t{d[subseq2]}\t{d[seqnum]}\t{d[hist400]}\t{res:.3f}\t{d[avlen]}\n'
                    .format(d=d, res=(float(d['res']) * float(d['hist400'])) / 63240))
print('done')
