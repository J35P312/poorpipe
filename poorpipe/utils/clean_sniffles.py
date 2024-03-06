import sys

for line in open(sys.argv[1]):
	if line[0] == "#":
		print(line.strip())
		continue
	
	if "SVTYPE=INS" in line or "SVTYPE=BND" in line or "SVTYPE:DUP/INS" in line:
		var=line.strip().split(";END=")
		before=var[0]
		after=var[-1].lstrip("0123456789")		
		line=before+after

	if "SVTYPE=BND" in line:
		var=line.strip().split(";SVLEN=1")
		before=var[0]
		after=var[-1]		
		line=before+after



	print(line.strip())
