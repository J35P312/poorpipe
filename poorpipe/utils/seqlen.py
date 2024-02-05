import sys
import math
import gzip

reads=500000

coverages={}
bin_list=list( range(0, 100, 1) )
for b in bin_list:
	coverages[b]=0

count=0
total=0
for line in gzip.open(sys.argv[1]):
	count+=1
	if total > reads:
		break

	if count == 2:
		rl=len(line)-1
		if rl > 100000:
			continue
		b=math.floor(rl/1000.0)

		total+=1
		coverages[b]+=1

	elif count == 4:
		count=0

print("# id: Read Length distribution")
print("# section_name: \'Read Length distribution\'")
print("# description: \'This output is described in the file header. Any MultiQC installation will understand it without prior configuration.\'")
print("# format: \'tsv\'")
print("# plot_type: \'linegraph\'")
print("# pconfig:")
print("#    id: 'custom_bargraph_w_header'")
print("#    title: Read Length distribution")
print("#    ylab: 'Fraction of reads'")
print("#    xlab: 'Read Lenght(kb)'")
print("#    height: 400")


for c in sorted(coverages.keys()):
	print("{}\t{}".format(c,coverages[c]/total))
