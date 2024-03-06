def align(input_folder,sample,family,ref,samtools,readlength_script,generate_ped_script):
	prefix=sample
	bam2fastq=[]

	fastq_list=open(f"{prefix}/parallel_fastq.txt","w")
	for line in open(f"{sample}/bam_input.txt"):
		b=line.split("/")[-1].split(".")[0]
		f=line.strip()
		bam2fastq.append(f"{samtools} fastq -T MM,ML {f} | gzip -c > {prefix}/fastq/{b}.fastq.gz" )

	bam2fastq="\n".join(bam2fastq)
	fastq_list.write(bam2fastq)

	cmd=f"""

parallel --tmpdir {sample}/tmp -j 20 < {prefix}/parallel_fastq.txt
cat {sample}/fastq/* > {prefix}/{prefix}.fastq.gz
python {readlength_script} {prefix}/{prefix}.fastq.gz > {prefix}/{prefix}_mqc.txt
minimap2 -R "@RG\\tID:{prefix}\\tSM:{prefix}" -a -y -t 16 --MD -x map-ont {ref} {prefix}/{prefix}.fastq.gz | {samtools} sort -T {sample}/tmp/{sample}.samtools -m 5G -@8 - > {prefix}/{prefix}.bam
{samtools} index -@16  {prefix}/{prefix}.bam
{samtools} idxstats {prefix}/{prefix}.bam > {prefix}/{prefix}.bam.bai.stats.txt
{samtools} stats {prefix}/{prefix}.bam > {prefix}/{prefix}.bam.stats.txt
python {generate_ped_script} {sample} {family} {prefix}/{prefix}.bam.bai.stats.txt > {prefix}/{prefix}.gender.txt
"""
	return(cmd)

