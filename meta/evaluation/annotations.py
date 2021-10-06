with open('ena_data_20210917-1328.txt') as f, open('annotations.csv', 'w') as a, open('readlengths.csv', 'w') as l:
    for line in f:
        if line.startswith('AC'):
            accession = line[2:].strip()[:-1]
        elif line.startswith('FT   gene'):
            span = line[9:].strip()
            if span.startswith('complement('):
                span = span[11:-1]
                strand = '-'
            else:
                strand = '+'
            start, end = span.split('..')
            print(accession, start, end, strand, sep=',', file=a)
        elif line.startswith('SQ'):
            parts = line[2:].strip().split(' ')
            assert parts[0] == "Sequence"
            print(accession, parts[1], sep=',', file=l)
            assert parts[2].startswith("BP")
