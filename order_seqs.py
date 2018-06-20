with open('subseqs.txt') as f:
    lines = list((str(len(lst[0])), lst[1], lst[2], lst[0]) for lst in (line.rstrip('\n').split('\t') for line in f.readlines()))
lines.sort(key=lambda x: int(x[1]), reverse=True)
with open('subseqs_ordered.txt', 'w') as f:
    for line in lines:
        f.write('\t'.join(line) + '\n')
print('done')
