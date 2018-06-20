with open('subseqs.txt') as f:
    lines = list(line.split('\t') for line in f.readlines())
lines.sort(key=lambda x: (len(x[0]), int(x[1])), reverse=True)
with open('subseqs_ordered.txt', 'w') as f:
    for line in lines:
        f.write('\t'.join(line))
print('done')
