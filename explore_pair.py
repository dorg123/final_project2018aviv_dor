import read
import sequence_funcs as sf
r = read.FastaReader('eztaxon_qiime_full.fasta')
sequences = list(r.data.values())
subseq1 = 'TTAGATACCCTGGTAGTC'
subseq2 = 'CCTACGGGAGGCAGCAGT'
lst1 = list(seq.find(subseq1) for seq in sequences)
lst2 = list(seq.find(subseq2) for seq in sequences)
non_mens1 = lst1.count(-1)
non_mens2 = lst2.count(-1)
mens1 = len(lst1) - non_mens1
mens2 = len(lst2) - non_mens2
avspot1 = (sum(lst1) + non_mens1) / mens1
avspot2 = (sum(lst2) + non_mens2) / mens2
while -1 in lst1:
    lst1.remove(-1)
while -1 in lst2:
    lst2.remove(-1)
lstz = list(abs(l1 - l2) if l1 != -1 and l2 != -1 else -1 for l1, l2 in zip(lst1, lst2))
while -1 in lstz:
    lstz.remove(-1)
print('Subseq1:')
print('Av. loc.: {:.2f}'.format(avspot1))
print('Max loc.: {}'.format(max(lst1)))
print('Min loc.: {}'.format(min(lst1)))
print(sf.number_histogram(lst1, 100))
print('Subseq2:')
print('Av. loc.: {:.2f}'.format(avspot2))
print('Max loc.: {}'.format(max(lst2)))
print('Min loc.: {}'.format(min(lst2)))
print(sf.number_histogram(lst2, 100))
print('diff:')
print('Av. len: {:.2f}'.format(sum(lstz) / len(lstz)))
print('Seq. num: {}'.format(len(lstz)))
print('Max len.: {}'.format(max(lstz)))
print('Min len.: {}'.format(min(lstz)))
print(sf.number_histogram(lstz, 100))
print('done')

