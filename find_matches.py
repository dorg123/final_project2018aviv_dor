#############################################
# W   W   A   RRR   N   N IIIII N   N  GGG  #
# W   W  A A  R  R  NN  N   I   NN  N G     #
# W W W AAAAA RRR   N N N   I   N N N G  GG #
#  W W  A   A R  R  N  NN   I   N  NN G   G #
#  W W  A   A R   R N   N IIIII N   N  GGG  #
#############################################
# This execution takes approximately 7-8 hours to complete.
import fasta
from time import time
# open the database, pairs and file
r = fasta.TabReader('pairs_by_hist400.txt')
db = fasta.FastaReader('eztaxon_qiime_full.fasta')
with open('pairs_resolution.txt', 'w') as f:
    f.write('subseq1\tsubseq2\tseqnum\tres\tavlen\thist400\thist500\n')
    for record in r.data:  # check for each pair
        seqs = []
        start = time()
        for seq in db.data.values():  # check for each sequence in the database
            if record['subseq1'] not in seq or record['subseq2'] not in seq:
                continue  # if the pair is not in the sequence, skip
            # if it is, add the sequence between the primers to the list
            subseq1loc, subseq2loc = tuple(sorted([(seq.find(record['subseq1']), len(record['subseq1'])),
                                                   (seq.find(record['subseq2']), len(record['subseq2']))]))
            seqs.append(seq[subseq1loc[0]:subseq2loc[0] + subseq2loc[1]])
        hist = []
        for seq in seqs:  # remove all the ones that appear more than once
            tmp = seq
            seqs.remove(tmp)
            if tmp in seqs:
                while tmp in seqs:
                    seqs.remove(tmp)
            else:
                hist.append(len(tmp))
        resolution = (len(hist) / len(seqs)) * 100  # and calculate the resolution
        f.write('{d[subseq1]}\t{d[subseq2]}\t{d[seqnum]}\t{res:.3f}\t{av:.1f}\t{d[hist400]}\t{d[hist500]}\n'
                .format(d=record, res=resolution, av=sum(hist) / len(hist)))
        end = time()
        print((end - start) * 1000)
print('done')
