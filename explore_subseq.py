import read
r = read.FastaReader('eztaxon_qiime_full.fasta')
subseq = 'TTAGATACCCTGGTAGTC'
lst = list(seq.find(subseq) for seq in r.data.values())
non_mens = lst.count(-1)
mens = len(lst) - non_mens
avspot = (sum(lst) + non_mens) / mens
while -1 in lst:
    lst.remove(-1)
print('Av. Loc.: {:.2f}'.format(avspot))
print('Max Loc.: {} ({})'.format(max(lst), len(max(r.data.values(), key=lambda x: len(x)))))
print('Min Loc.: {} ({})'.format(min(lst), len(min(r.data.values(), key=lambda x: len(x)))))
print('done')

