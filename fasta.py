class FileReader:
    def __init__(self, filename, accession=None):
        self._filename = filename
        self._read(accession)

    @property
    def filename(self):
        return self._filename

    @property
    def data(self):
        return self._data

    def __getitem__(self, item):
        try:
            return self._data[item]
        except KeyError:
            return None
        except IndexError:
            return None

    def __str__(self):
        return 'read.FileReader: {}, {} entries' \
            .format(self._filename, len(self._data))

    def __repr__(self):
        return 'read.FileReader({})<length: {}>' \
            .format(self._filename, len(self._data))


class FastaReader(FileReader):
    def _read(self, accession=None, read_head=lambda x: x):
        with open(self._filename, 'r') as f:
            self._data = fasta_reader(f, accession, read_head)

    @staticmethod
    def read_head(head, raw=False):
        parts = head.split('|')
        if raw:
            return dict(zip(range(len(parts)), parts))

        def find_all(s, sub):
            index = [s.find(sub)]
            while index[-1] != -1:
                index.append(s.find(sub, index[-1] + 1))
            index.pop()
            return index

        last = parts[-1]
        d = dict((parts[i], parts[i + 1]) for i in range(0, len(parts) - 1, 2))
        if '=' in last:
            separators = [0] + list(last.rfind(' ', 0, i) for i in find_all(last, '=')) + [len(last)]
            parts = list(last[separators[i - 1]:separators[i]].strip() for i in range(1, len(separators)))
            d.update(dict((part[:part.find('=')], part.split('=')[1]) if '=' in part else (0, part) for part in parts))
        return d

    def __str__(self):
        return 'fasta.FastaReader: {}, {} entries' \
            .format(self.filename, len(self.data))

    def __repr__(self):
        return 'fasta.FastaReader({})<length: {}>' \
            .format(self.filename, len(self.data))


class TabReader(FileReader):
    @property
    def head(self):
        return self._head

    def _read(self, accession=None):
        with open(self._filename, 'r') as f:
            lines = f.readlines()
        self._head, self._data = tab_reader(lines, True)

    def __str__(self):
        return 'read.TabReader: {} (head: {}), {} entries' \
            .format(self.filename, ', '.join(self.head), len(self.data))

    def __repr__(self):
        return 'read.TabReader({})<length: {}; head: {}>' \
            .format(self.filename, len(self.data), self.head)


class MapReader(FileReader):
    def _read(self, accession=None):
        with open(self._filename, 'r') as f:
            lines = f.readlines()
        self._data = tab_reader(lines, False)

    def __str__(self):
        return 'read.MapReader: {}, {} entries' \
            .format(self.filename, len(self.data))

    def __repr__(self):
        return 'read.MapReader({})<length: {}>' \
            .format(self.filename, len(self.data))


def fasta_reader(lines, accession=None, read_head=lambda x: x):
    d = dict()
    info = None
    flag = False
    for line in lines:
        if line.startswith('>'):
            if flag:
                break
            info = read_head(line.rstrip('\n').lstrip('>'))
            if info == accession:
                flag = True
        else:
            if accession is None or flag:
                d[info] = d.get(info, '') + line.rstrip('\n')
    return d


def tab_reader(lines, common_head=False):
    it = iter(lines)
    if common_head:
        head = next(it).rstrip('\n').split('\t')
        return head, list(dict(zip(head, line.rstrip('\n').split())) for line in it)
    return list(line.rstrip('\n').split() for line in it)


def fastq_reader(lines):
    d = dict()
    head = ''
    reading = False
    for line in lines:
        if line.startswith('@'):
            head = line.lstrip('@')
            reading = True
            d[head] = '', ''
        elif line.startswith('+'):
            reading = False
        else:
            s, a = d[head]
            if reading:
                s += line
            else:
                a += line
            d[head] = s, a
    return d


def plain_histogram(lst):
    return dict((key, lst.count(key)) for key in sorted(set(lst)))


def number_histogram(lst, size):
    return plain_histogram(list(map(lambda x: round(x / size) * size, lst)))


def subsequences(sequence, length):
    if len(sequence) <= length:
        return [sequence] if len(sequence) == length else []
    else:
        for i in range(len(sequence) - length):
            yield sequence[i:i + length]


def complementary(dna_seq):
    switch_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(switch_dict[b] for b in dna_seq)
