#############################################
# W   W   A   RRR   N   N IIIII N   N  GGG  #
# W   W  A A  R  R  NN  N   I   NN  N G     #
# W W W AAAAA RRR   N N N   I   N N N G  GG #
#  W W  A   A R  R  N  NN   I   N  NN G   G #
#  W W  A   A R   R N   N IIIII N   N  GGG  #
#############################################
# This execution takes approximately 30 minutes to complete.
import fasta


# explore_pair: takes two subsequences and a sequence database and checks for average length, hit number and histogram
def explore_pair(subseq1, subseq2, sequences):
    lst1 = [seq.find(subseq1) for seq in sequences]  # find subseq1 in the database
    lst2 = [seq.find(subseq2) for seq in sequences]  # find subseq2 in the database
    lst1 = [i for i in lst1 if i != -1]  # filter the -1s
    lst2 = [i for i in lst2 if i != -1]
    lstz = [abs(l1 - l2) if l1 != -1 and l2 != -1 else -1 for l1, l2 in zip(lst1, lst2)]  # calculate the distance
    lstz = [i for i in lstz if i != -1]  # filter the -1s
    return sum(lstz) / len(lstz), len(lstz), fasta.number_histogram(lstz, 100)  # return results


# read the database and filter for 40000 hits at least
reader = fasta.TabReader('subseqs.txt')
filtered = [seq for seq in reader.data if int(seq['seqnum']) > 40000]
# create pairs
pairs = []
for subseq1 in filtered:
    for subseq2 in filtered:
        if 600 > float(subseq1['avloc']) - float(subseq2['avloc']) > 400:  # filter by length between 400 and 600
            pairs.append((subseq1['subseq'], subseq2['subseq']))
# go over the database with each pair
reader = fasta.FastaReader('eztaxon_qiime_full.fasta')
sequences = reader.data.values()
with open('pairs.txt', 'w') as f:
    f.write('subseq1\tsubseq2\tseqnum\tavlen\thist400\thist500\n')
    for subseq1, subseq2 in pairs:
        av, seqnum, hist = explore_pair(subseq1, subseq2, sequences)  # explore_pair with the pairs
        f.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(subseq1, subseq2, seqnum, av, hist[400], hist[500]))
        print('*')
print('done')
