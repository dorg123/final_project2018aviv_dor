#############################################
# W   W   A   RRR   N   N IIIII N   N  GGG  #
# W   W  A A  R  R  NN  N   I   NN  N G     #
# W W W AAAAA RRR   N N N   I   N N N G  GG #
#  W W  A   A R  R  N  NN   I   N  NN G   G #
#  W W  A   A R   R N   N IIIII N   N  GGG  #
#############################################
# This execution takes approximately 1 hour to complete.
import fasta
# read database
r = fasta.FastaReader('eztaxon_qiime_full.fasta')
seq = r['130648']  # e.coli 16s
# create subsequences
subseqs = list(fasta.subsequences(seq, 18)) + list(fasta.subsequences(seq, 19)) + list(fasta.subsequences(seq, 20))\
          + list(fasta.subsequences(seq, 21)) + list(fasta.subsequences(seq, 22))
# iterate over the subsequences
lines = []
for subseq in subseqs:
    # for each one, iterate over the database and seek their position in each sequence
    lst = list(seq.find(subseq) for seq in r.data.values())
    # count how many sequences contain the specific subsequence
    mens = len(lst) - lst.count(-1)
    # calculate the average location
    avspot = (sum(lst) + lst.count(-1)) / mens
    # and save the data
    lines.append((str(len(subseq)), str(mens), str(avspot), subseq))
    print('*')
# sort the data
lines.sort(key=lambda x: int(x[1]), reverse=True)
with open('subseqs.txt', 'w') as f:  # and write it to a file
    f.write('len\tseqnum\tavloc\tsubseq\n')
    for line in lines:
        f.write('\t'.join(line) + '\n')
print('done')

