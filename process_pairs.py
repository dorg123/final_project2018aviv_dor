import read
reader = read.TabReader('pairs.txt')
data = [d for d in sorted(reader.data, key=lambda x: int(x['SeqNum']), reverse=True) if int(d['Hist400']) > 30000]
with open('pairs_by_seqnum.txt', 'w') as f:
    f.write('\t'.join(reader.head) + '\n')
    for d in data:
        f.write('{}\t{}\t{}\t{:.1f}\t{}\t{}\n'.format(d['SubSeq1'], d['SubSeq2'], d['SeqNum'],
                                                      float(d['AvLoc']), d['Hist400'], d['Hist500']))
print('done')
