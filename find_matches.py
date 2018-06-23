#############################################
# W   W   A   RRR   N   N IIIII N   N  GGG  #
# W   W  A A  R  R  NN  N   I   NN  N G     #
# W W W AAAAA RRR   N N N   I   N N N G  GG #
#  W W  A   A R  R  N  NN   I   N  NN G   G #
#  W W  A   A R   R N   N IIIII N   N  GGG  #
#############################################
# This execution takes approximately 5-10 minutes to complete.
import read
from time import time
r = read.TabReader('pairs_by_hist400.txt')
db = read.FastaReader('eztaxon_qiime_full.fasta')
with open('pairs_resolution.txt', 'w') as f:
    f.write('SubSeq1\tSubSeq2\tSeqNum\tRes\tAvLen\tHist400\tHist500\n')
    for record in r.data:
        seqs = []
        start = time()
        for seq in db.data.values():
            if record['SubSeq1'] not in seq or record['SubSeq2'] not in seq:
                continue
            subseq1loc, subseq2loc = tuple(sorted([(seq.find(record['SubSeq1']), len(record['SubSeq1'])),
                                                   (seq.find(record['SubSeq2']), len(record['SubSeq2']))]))
            seqs.append(seq[subseq1loc[0]:subseq2loc[0] + subseq2loc[1]])
        hist = set(seqs)
        resolution = (len(hist) / len(seqs)) * 100
        f.write('{d[SubSeq1]}\t{d[SubSeq2]}\t{d[SeqNum]}\t{res:.2f}\t{d[AvLoc]}\t{d[Hist400]}\t{d[Hist500]}\n'
                .format(d=record, res=resolution))
        end = time()
        print((end - start) * 1000)
print('done')
