def cnvpytor_call(bam,sample,ref,pytor):
	cmd=f"""{pytor} -root {sample}/{sample}.pytor -rd {bam}
{pytor} -root {sample}/{sample}.pytor -gc {ref}
{pytor} -root {sample}/{sample}.pytor -his 2000 200000
{pytor} -root {sample}/{sample}.pytor -partition 2000 200000
{pytor} -root {sample}/{sample}.pytor -call 2000 > {sample}/{sample}.pytor.out
{pytor} -root {sample}/{sample}.pytor -call 200000 > {sample}/{sample}.pytor.200k.out

{pytor} -root {sample}/{sample}.pytor -view 2000 <<ENDL
set print_filename {sample}/{sample}.pytor.vcf
print calls
ENDL

"""

	return(cmd)

def cnvpytor_baf(bam,sample,ref,pytor):


	manhattan=f"""
{pytor} -root {sample}/{sample}.pytor -view 200000 <<ENDL
set style bmh
set rd_use_mask
set file_titles {sample}
manhattan
save {sample}/pytor/{sample}.manhattan.png
ENDL
"""

	chr_plots=[]
	chromosomes=list(range(1,23))+["X","Y"]

	for chr in chromosomes:
		chr_plots.append(f"""
{pytor} -root {sample}/{sample}.pytor -view 200000 <<ENDL
set style classic
set rd_use_mask
set file_titles {sample}
set panels rd likelihood snp
set markersize 0.2
{chr}

save {sample}/pytor/{sample}.{chr}.png
ENDL
""")
	chr_plots="\n".join(chr_plots)


	cmd=f"""
{pytor} -root {sample}/{sample}.pytor -snp {sample}/{sample}.pytor.vcf.gz
{pytor} -root {sample}/{sample}.pytor -view 200000 <<ENDL
{pytor} -root file.pytor -mask_snps
{pytor} -root file.pytor -baf 200000

{chr_plots}

{manhattan}

	"""
	return(cmd)

