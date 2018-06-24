#############################################
# W   W   A   RRR   N   N IIIII N   N  GGG  #
# W   W  A A  R  R  NN  N   I   NN  N G     #
# W W W AAAAA RRR   N N N   I   N N N G  GG #
#  W W  A   A R  R  N  NN   I   N  NN G   G #
#  W W  A   A R   R N   N IIIII N   N  GGG  #
#############################################
# This execution takes approximately 30 minutes to complete.
import fasta


def explore_pair(subseq1, subseq2, sequences):
    lst1 = [seq.find(subseq1) for seq in sequences]
    lst2 = [seq.find(subseq2) for seq in sequences]
    lst1 = [i for i in lst1 if i != -1]
    lst2 = [i for i in lst2 if i != -1]
    lstz = [abs(l1 - l2) if l1 != -1 and l2 != -1 else -1 for l1, l2 in zip(lst1, lst2)]
    lstz = [i for i in lstz if i != -1]
    return sum(lstz) / len(lstz), len(lstz), fasta.number_histogram(lstz, 100)


reader = fasta.TabReader('subseqs.txt')
filtered = [seq for seq in reader.data if int(seq['seqnum']) > 40000]
pairs = []
for subseq1 in filtered:
    for subseq2 in filtered:
        if 600 > float(subseq1['avloc']) - float(subseq2['avloc']) > 400:
            pairs.append((subseq1['subseq'], subseq2['subseq']))
reader = fasta.FastaReader('eztaxon_qiime_full.fasta')
sequences = reader.data.values()
with open('pairs.txt', 'w') as f:
    f.write('subseq1\tsubseq2\tseqnum\tavlen\thist400\thist500\n')
    for subseq1, subseq2 in pairs:
        av, seqnum, hist = explore_pair(subseq1, subseq2, sequences)
        f.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(subseq1, subseq2, seqnum, av, hist[400], hist[500]))
        print('*')
print('done')
