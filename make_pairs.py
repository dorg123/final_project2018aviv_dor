import read
import sequence_funcs as sf


def explore_pair(subseq1, subseq2, sequences):
    lst1 = [seq.find(subseq1) for seq in sequences]
    lst2 = [seq.find(subseq2) for seq in sequences]
    lst1 = [i for i in lst1 if i != -1]
    lst2 = [i for i in lst2 if i != -1]
    lstz = [abs(l1 - l2) if l1 != -1 and l2 != -1 else -1 for l1, l2 in zip(lst1, lst2)]
    lstz = [i for i in lstz if i != -1]
    return sum(lstz) / len(lstz), len(lstz), sf.number_histogram(lstz, 100)


reader = read.MapReader('subseqs.txt')
filtered = [seq for seq in reader.data if int(seq[1]) > 40000]
pairs = []
for subseq1 in filtered:
    for subseq2 in filtered:
        if 600 > float(subseq1[2]) - float(subseq2[2]) > 400:
            pairs.append((subseq1[3], subseq2[3]))
reader = read.FastaReader('eztaxon_qiime_full.fasta')
sequences = reader.data.values()
with open('pairs.txt', 'w') as f:
    f.write('SubSeq1\tSubSeq2\tSeqNum\tAvLoc\tHist400\tHist500\n')
    for subseq1, subseq2 in pairs:
        av, seqnum, hist = explore_pair(subseq1, subseq2, sequences)
        f.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(subseq1, subseq2, seqnum, av, hist[400], hist[500]))
print('done')
