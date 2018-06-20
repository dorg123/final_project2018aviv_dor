import read
import sequence_funcs as sf
import time as tm

prev = tm.time()


def time():
    global prev
    tmp = prev
    prev = tm.time()
    return 1000 * (tm.time() - tmp)


r = read.FastaReader('eztaxon_qiime_full.fasta')
seq = r['130648']
time()
subseqs = list(sf.subsequences(seq, 18)) + list(sf.subsequences(seq, 19)) + list(sf.subsequences(seq, 20))\
          + list(sf.subsequences(seq, 21)) + list(sf.subsequences(seq, 22))
print(time())
print(len(subseqs))
d = dict()
for subseq in subseqs:
    time()
    d[subseq] = list(seq.find(subseq) for seq in r.data.values())
    print(time())
with open('subseqs.txt', 'w') as f:
    for subseq, lst in d.items():
        mens = len(lst) - lst.count(-1)
        avspot = (sum(lst) + lst.count(-1)) / mens
        f.write('{}\t{}\t{:.2f}\n'.format(subseq, mens, avspot))
print('done')

