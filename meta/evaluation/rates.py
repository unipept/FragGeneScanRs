#!/usr/bin/env python
from collections import namedtuple
from operator import itemgetter

Annotation = namedtuple('Annotation', ['name', 'start', 'end', 'strand'])
Event = namedtuple('Event', ['location', 'annotation', 'prediction'])

def parse_readlengths(filename):
	with open(filename) as f:
		for line in f:
			name, length = line.strip().split(',')
			yield name, int(length)

def parse_csv(filename):
	with open(filename) as f:
		for line in f:
			name, start, end, strand = line.strip().split(',')
			yield Annotation(name, int(start), int(end), 1 if strand == '+' else -1)

def events(read_lengths, annotations, predictions):
	names = dict()
	total = 0
	for name, length in parse_readlengths(read_lengths):
		names[name] = total
		total += length

	for name, start, end, strand in parse_csv(annotations):
		yield Event(names[name] + start, strand, None)
		yield Event(names[name] + end, 0, None)

	for name, start, end, strand in parse_csv(predictions):
		yield Event(names[name] + start, None, strand)
		yield Event(names[name] + end, None, 0)

def rates(read_lengths, annotations, predictions):
	rates = dict(tp=0, fp=0, tn=0, fn=0)
	cl = 0 # current location
	ca = 0 # current annotation
	cp = 0 # current prediction
	for l, a, p in sorted(events(read_lengths, annotations, predictions), key=itemgetter(0)):
		if ca == 0 and cp == 0:
			rates['tn'] += l - cl
		elif ca == cp:
			rates['tp'] += l - cl
		elif ca == 0 and cp != 0:
			rates['fp'] += l - cl
		elif ca != 0 and cp == 0:
			rates['fn'] += l - cl
		else: # different strands
			rates['fp'] += l - cl

		if a is not None: ca = a
		if p is not None: cp = p
		cl = l
	return rates

row = '{:<15} {:>10} {:>10} {:>10} {:>10} {:>10.4}{:>10.4}{:>10.4}{:>10.4}{:>10.4}'
print(row.format('tool', 'TP', 'FP', 'TN', 'FN', 'precision', 'sensitivity', 'specificity', 'NPV', 'MCC'))
for tool in ['FGS', 'FGS+', 'prodigal', 'FGSrs']:
	r = rates('readlengths.csv', 'annotations.csv', f'{tool}.csv')
	tp, fp, tn, fn = r['tp'], r['fp'], r['tn'], r['fn']
	print(row.format(tool, tp, fp, tn, fn, 100 * tp / (tp + fp), 100 * tp / (tp + fn), 100 * tn / (tn + fp), 100 * tn / (tn + fn), 100 * (tp * tn - fp * fn) / ((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn))**0.5))
