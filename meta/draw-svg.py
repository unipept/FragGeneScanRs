#!/usr/bin/env python
import sys
import drawSvg as draw # pip install drawSvg

class Line:
	def __init__(self, tool, typ, reads, threads, time, memory, resmem):
		self.tool = tool
		self.typ = typ
		self.reads = int(reads)
		self.threads = int(threads)
		self.time = int(time) / 1000
		self.memory = int(memory) / 1024 / 1024
		self.resmem = int(resmem) / 1024 / 1024

values = []
with open('benchmark.csv') as csv:
	_ = next(csv) # header
	for line in csv:
		values.append(Line(*line.split(',')))

threadnum = 16
radius = 1

color = { 'FGS': '#4e79a7', 'FGS+': '#f28e2b', 'FGSrs': '#e15759', 'FGSrsu': '#76b7b2' }

width = 100
height = 100
parallel = draw.Drawing(width, height)
parallel.append(draw.Rectangle(0, 0, width, height, fill='white'))
parallel.append(draw.Line(0, 0, width, height, stroke='grey', stroke_dasharray="2, 2"))

bases = {}
for v in values:
	if v.typ != 'long': continue
	if v.threads == 1:
		bases[v.tool] = v.time
	parallel.append(draw.Circle(v.threads * width / threadnum,
	                bases[v.tool] * height / threadnum / v.time,
	                radius, fill=color[v.tool]))

parallel.saveSvg('parallel-efficiency.svg')

scale = 10
width = 100
height = max(v.memory for v in values) / scale
memory_usage = draw.Drawing(width, height)
memory_usage.append(draw.Rectangle(0, 0, width, height, fill='white'))

for v in values:
	if v.typ != 'long': continue
	memory_usage.append(draw.Circle(v.threads * width / threadnum, v.memory / scale, radius, fill=color[v.tool]))

memory_usage.saveSvg('memory-usage.svg')

print("""
| Short reads      |  1 thread | 2 threads | 4 threads | 8 threads | 16 threads |
|:-----------------|----------:|----------:|----------:|----------:|-----------:|
| FragGeneScan     | {fgs_[1]} r/s | {fgs_[2]} r/s | {fgs_[4]} r/s | {fgs_[8]} r/s | {fgs_[16]} r/s | 
| FragGeneScanPlus | {fgsp[1]} r/s | {fgsp[2]} r/s | {fgsp[4]} r/s | {fgsp[8]} r/s | / | 
| FragGeneScanRs   | {fgsr[1]} r/s | {fgsr[2]} r/s | {fgsr[4]} r/s | {fgsr[8]} r/s | {fgsr[16]} r/s | 
""".format(fgs_={ v.threads: round(v.reads / v.time) for v in values if v.typ == "short" and v.tool == "FGS" },
           fgsp={ v.threads: round(v.reads / v.time) for v in values if v.typ == "short" and v.tool == "FGS+" },
           fgsr={ v.threads: round(v.reads / v.time) for v in values if v.typ == "short" and v.tool == "FGSrs" }))

print("""
| Long reads       |  1 thread | 2 threads | 4 threads | 8 threads | 16 threads |
|:-----------------|----------:|----------:|----------:|----------:|-----------:|
| FragGeneScan     | {fgs_[1]} r/s | {fgs_[2]} r/s | {fgs_[4]} r/s | {fgs_[8]} r/s | {fgs_[16]} r/s | 
| FragGeneScanPlus | {fgsp[1]} r/s | {fgsp[2]} r/s | {fgsp[4]} r/s | {fgsp[8]} r/s | / | 
| FragGeneScanRs   | {fgsr[1]} r/s | {fgsr[2]} r/s | {fgsr[4]} r/s | {fgsr[8]} r/s | {fgsr[16]} r/s | 
""".format(fgs_={ v.threads: round(v.reads / v.time) for v in values if v.typ == "long" and v.tool == "FGS" },
           fgsp={ v.threads: round(v.reads / v.time) for v in values if v.typ == "long" and v.tool == "FGS+" },
           fgsr={ v.threads: round(v.reads / v.time) for v in values if v.typ == "long" and v.tool == "FGSrs" }))

print("""
| Complete genome  |  1 thread |
|:-----------------|----------:|
| FragGeneScan     | {fgs_} s |
| FragGeneScanPlus | {fgsp} s |
| FragGeneScanRs   | {fgsr} s |
""".format(fgs_= next((round(v.time, 3) for v in values if v.typ == "complete" and v.tool == "FGS")),
           fgsp= next((round(v.time, 3) for v in values if v.typ == "complete" and v.tool == "FGS+")),
           fgsr= next((round(v.time, 3) for v in values if v.typ == "complete" and v.tool == "FGSrs"))))

print("""
| Short reads      |  1 thread | 2 threads | 4 threads | 8 threads |
|:-----------------|----------:|----------:|----------:|----------:|
| FragGeneScanPlus | {fgsp[1]} r/s | {fgsp[2]} r/s | {fgsp[4]} r/s | {fgsp[8]} r/s |
| FragGeneScanRs   | {fgsr[1]} r/s | {fgsr[2]} r/s | {fgsr[4]} r/s | {fgsr[8]} r/s |
""".format(fgsp={ v.threads: round(v.reads / v.time) for v in values if v.tool == "stdout FGS+" },
           fgsr={ v.threads: round(v.reads / v.time) for v in values if v.tool == "stdout FGSrsu" }))

def draw_absolute(values, length):
	fgs = next((round(v.reads / v.time) for v in values if v.typ == length and v.tool == "FGS"))
	fgsp = next((round(v.reads / v.time) for v in values if v.typ == length and v.tool == "FGS+"))
	fgsrs = next((round(v.reads / v.time) for v in values if v.typ == length and v.tool == "FGSrs"))
	reads = next((v.reads for v in values if v.typ == length))
	width = 200
	height = 32
	absolute = draw.Drawing(width, height)
	absolute.append(draw.Rectangle(0, 0, width, height, fill='white'))
	maximum = max([fgs, fgsp, fgsrs])
	zeros = 10 ** (len(str(maximum)) - 2)
	maximum = maximum // zeros * zeros + zeros
	absolute.append(draw.Rectangle(0, 0, fgsrs * width / maximum, 10, fill=color['FGSrs']))
	absolute.append(draw.Rectangle(0, 11, fgsp * width / maximum, 10, fill=color['FGS+']))
	absolute.append(draw.Rectangle(0, 22, fgs * width / maximum, 10, fill=color['FGS']))
	absolute.saveSvg(f'absolute-{length}.svg')
	print(f'absolute {length} maximum: {maximum}')

draw_absolute(values, "short")
draw_absolute(values, "long")
