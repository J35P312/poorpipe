import sys

sample=sys.argv[1]
family=sys.argv[2]
idx_stat=sys.argv[3]

chromosomes={}

X=False
Y=False
for line in open(idx_stat):
	content=line.strip().split()
	chromosomes[content[0]]=float(content[2])
	if content[0] == "chrY" or content[0] == "Y":
		Y=content[0]
	if content[0] == "chrX" or content[0] == "X":
		X=content[0]

sex=2
if X and Y:
	if chromosomes[X]:
		if chromosomes[Y]/chromosomes[X] > 0.2:
			sex=1

if sex == 2:
	print("female")
else:
	print("male")


ped=open(f"{sample}/{family}.fam","w")
ped.write(f"#family_id\tsample_id\tfather\tmother\tsex\tphenotype\n")
ped.write(f"{family}\t{sample}\t0\t0\t{sex}\t2")
ped.close()

