#!/usr/bin/env python
import pysam # pip install pysam
from operator import attrgetter

class DataClass:
	__slots__ = ()

	def __init__(self, *args):
		for name, value in zip(self.__slots__, args):
			setattr(self, name, value)

	def __repr__(self):
		args = ', '.join(f'{repr(getattr(self, name))}' for name in self.__slots__)
		return f'{type(self).__name__}({args})'

	def aslist(self):
		return [type(self)] + [getattr(self, name) for name in self.__slots__]

	def __eq__(self, other):
		return self.aslist() == other.aslist()

	def __lt__(self, other):
		return self.aslist() < other.aslist()


# ============================================================
# Reading genes and sequence lengths from annotations

class Gene(DataClass):
	__slots__ = ('accession', 'start', 'end', 'strand', 'translation')

genes = []
sequences = {}

with open('ena_data_20210917-1328.txt') as f:
	in_translation = False
	for line in f:
		line = line.strip()
		if not in_translation:
			if line.startswith('AC'):
				accession = line[2:].strip()[:-1]
			elif line.startswith('FT   gene'):
				span = line[9:].strip()
				if span.startswith('complement('):
					span = span[11:-1]
					strand = -1
				else:
					strand = 1
				start, end = span.split('..')
			elif line.startswith('FT                   /translation="'):
				in_translation = True
				translation = ""
			elif line.startswith('SQ'):
				parts = line[2:].strip().split(' ')
				assert parts[0] == "Sequence"
				sequences[accession] = int(parts[1])
				assert parts[2].startswith("BP")
		if in_translation:
			translation += line.split(' ')[-1]
			if translation.endswith('"'):
				_, translation, _ = translation.split('"')
				genes.append(Gene(accession, int(start), int(end), strand, translation))
				in_translation = False

def get_genes_in_region(name, start, end):
	for gene in genes:
		if gene.accession != name: continue
		if gene.end < start: continue
		if end < gene.start: continue
		yield Gene(name, max(start, gene.start), min(end, gene.end), gene.strand, gene.translation)

# ============================================================
# Reading predictions from FGS*

class Prediction(DataClass):
	__slots__ = ('origin', 'start', 'end', 'strand', 'sequence')

predictions = dict()

def to_strand(char):
	return 1 if char == '+' else -1

for predictor in ['FGS', 'FGS+', 'FGSrs']:
	with open(f'{predictor}2000.faa') as f:
		header = next(f)
		sequence = ""
		for line in f:
			if line.startswith('>'):
				origin, start, end, strand = header[1:].split('_')
				current = predictions.setdefault(predictor, {})
				preds = current.setdefault(origin, [])
				preds.append(Prediction(origin, int(start), int(end), to_strand(strand), sequence))

				header = line.strip()
				sequence = ""
			else:
				sequence += line.strip()

		origin, start, end, strand = header[1:].split('_')
		current = predictions.setdefault(predictor, {})
		preds = current.setdefault(origin, [])
		preds.append(Prediction(origin, int(start), int(end), to_strand(strand), sequence))

# ============================================================
# Reading simulated read offsets from SAM-file

class Offset(DataClass):
	__slots__ = ('reference', 'start', 'end', 'strand', 'name', 'segment')

offsets = {}

body = '{:<10}{:>8.2%}{:>8.2%}{:>8.2%}{:>8.2%}{:>8.2%}{:>8.2%}{:>8.2%}{:>8.2%}{:>8.2}'
head = '{:<10}{:>8.4s}{:>8.4s}{:>8.4s}{:>8.4s}{:>8.4s}{:>8.4s}{:>8.4s}{:>8.4s}{:>8.4s}'
print(head.format('tool', 'TP', 'FP', 'TN', 'FN', 'precision', 'sensitivity', 'specificity', 'NPV', 'MCC'))

for predictor in predictions:
	tp, fp, tn, fn, ds = 0, 0, 0, 0, 0
	with open('2000.sam') as f:
		samfile = pysam.AlignmentFile(f)
		for segment in samfile:
			origin = segment.get_tag('oR').split('|')[1]
			preds = predictions[predictor].get(segment.query_name, [])
			local_genes = list(get_genes_in_region(origin, int(segment.reference_start), int(segment.reference_end)))
			#print('GOOOOOOO', segment.reference_start, segment.reference_end, local_genes, preds)
			for (sim, ref) in segment.get_aligned_pairs():
				expected = { gene.strand for gene in local_genes if ref is not None and gene.start <= ref and ref < gene.end }
				predicted = { -prediction.strand if segment.is_reverse else prediction.strand for prediction in preds if sim is not None and prediction.start <= sim and sim < prediction.end }
				#print(ref, sim, expected, predicted, end=' ')
				if len(expected) == 0 and len(predicted) == 0:
					#print('tn')
					tn += 1
				elif len(expected) == 0 and len(predicted) != 0:
					#print('fp')
					fp += 1
				elif len(expected) != 0 and len(predicted) == 0:
					#print('fn')
					fn += 1
				elif set(expected) == set(predicted):
					#print('tp')
					tp += 1
				else:
					#print('ds')
					fp += 1
					ds += 1

	t = tp + fp + tn + fn
	print(body.format(
		predictor,
		tp / t,
		fp / t,
		tn / t,
		fn / t,
		tp / (tp + fp),
		tp / (tp + fn),
		tn / (tn + fp),
		tn / (tn + fn),
		(tp * tn - fp * fn) / ((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn))**0.5
	), f'{ds/t:.2%}')

