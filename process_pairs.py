import fasta
# read the pairs
reader = fasta.TabReader('pairs.txt')
data = [d for d in reader.data if int(d['hist400']) > 30000]  # filter by hit number greater than 30,000
with open('pairs_by_hist400.txt', 'w') as f:
    f.write('\t'.join(reader.head) + '\n')
    for d in data:
        f.write('{d[subseq1]}\t{d[subseq2]}\t{d[seqnum]}\t{avlen:.1f}\t{d[hist400]}\t{d[hist500]}\n'
                .format(d=d, avlen=float(d['avlen'])))  # and write the pairs
print('done')
