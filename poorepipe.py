from slurmpy import Slurm
import sys
import os
import glob

account="development"
samtools="singularity exec --bind /home/proj/ singularity/samtools_1.19--h50ea8bc_0.sif samtools"
sniffles="singularity exec --bind /home/proj/ singularity/sniffles_1.0.12--h8b12597_1.sif sniffles"
picard="singularity exec --bind /home/proj/ singularity/picard_3.1.1--hdfd78af_0.sif picard"
pytor="singularity exec --bind /home/proj singularity/cnvpytor_1.3.1--pyhdfd78af_1.sif cnvpytor"
deepvariant="singularity exec --bind /home/proj/ singularity/deepvariant_latest.sif"
bcftools="singularity exec --bind /home/proj/ singularity/bcftools_1.19--h8b25389_0.sif bcftools"
nanostats="singularity exec --bind /home/proj/ singularity/nanostat_1.6.0--pyhdfd78af_0.sif NanoStat"
fastqc="singularity exec --bind /home/proj/ singularity/fastqc_0.12.1--hdfd78af_0.sif fastqc"
multiqc="singularity exec --bind /home/proj/ singularity/multiqc_1.18--pyhdfd78af_0.sif multiqc"
whatshap="singularity exec --bind /home/proj/ singularity/whatshap_2.1--py39h1f90b4d_0.sif whatshap"

trgt="singularity exec --bind /home/proj/ singularity/trgt_0.4.0.sif trgt"
trgt_repeats="reeats_hg19.bed"

ref="/home/proj/production/rare-disease/references/references_12.0/grch37_homo_sapiens_-d5-.fasta"

#{samtools} merge -f {prefix}/{prefix}.bam {input_folder}/*bam
#{samtools} fastq -T MM,ML {prefix}/{prefix}.bam | gzip -c > {prefix}/{prefix}.fastq.gz
def align(input_folder,sample,ref,samtools):
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

parallel -j 20 < {prefix}/parallel_fastq.txt
cat {sample}/fastq/* > {prefix}/{prefix}.fastq.gz
minimap2 -R "@RG\\tID:{prefix}\\tSM:{prefix}" -a -t 16 --MD -x map-ont {ref} {prefix}/{prefix}.fastq.gz | {samtools} sort -m 5G -@8 - > {prefix}/{prefix}.bam
{samtools} index -@16  {prefix}/{prefix}.bam
{samtools} stats {prefix}/{prefix}.bam > {prefix}/{prefix}.bam.stats.txt
"""
	return(cmd)

def sniffles_sv(bam,sample,sniffles):
	cmd=f"""{sniffles} -m {bam} -v {sample}/{sample}.sniffles.vcf -l 100 -t 16 -s 3 --genotype --cluster"""
	return(cmd)

def target_type(bam,sample,ref,trgt_repeats,trgt):
	cmd=f"{trgt} --genome {ref} --repeats {trgt_repeats} --reads {bam} --output-prefix {sample}/{sample}.trgt"
	return(cmd)

def picard_qc(bam,sample,ref,picard):
	cmd=f"{picard}  -Xmx20G CollectWgsMetrics --INPUT {bam} --OUTPUT {sample}/{sample}.wgsmetrics.txt --REFERENCE_SEQUENCE {ref} --COUNT_UNPAIRED true --MINIMUM_BASE_QUALITY 1 --VALIDATION_STRINGENCY SILENT"
	return(cmd)

def cnvpytor_call(bam,sample,ref,pytor):
	cmd=f"""{pytor} -root {sample}/{sample}.pytor -rd {bam}
{pytor} -root {sample}/{sample}.pytor -gc {ref}
{pytor} -root {sample}/{sample}.pytor -his 2000 200000
{pytor} -root {sample}/{sample}.pytor -partition 2000 200000
{pytor} -root {sample}/{sample}.pytor -call 2000 > {sample}/{sample}.pytor.out"""

	return(cmd)

def deepvariant_call_chr(bam,sample,ref,deepvariant,chromosome):
	cmd=f"""
TMPDIR={sample}/tmp/
{deepvariant} run_deepvariant --output_vcf {sample}/{sample}/dv.{chromosome}.vcf --model_type ONT_R104 --num_shards 16 --reads {bam} --ref {ref} --regions {chromosome} --intermediate_results_dir {sample}/tmp/{chromosome}"""
	return(cmd)

def bcftools_concat(sample,bcftools):

	bcftools_cmd=f"""
{bcftools} concat {sample}/{sample}/dv.*vcf | {bcftools} sort -O z - --temp-dir {sample}/tmp/ > {sample}/{sample}.vcf.gz
{bcftools} index {sample}/{sample}.vcf.gz
{bcftools} stats {sample}/{sample}.vcf.gz > {sample}/{sample}.txt
	
"""
	return(bcftools_cmd)

def whatshap_phase(bam,sample,ref,bcftools,whatshap):
	cmd=f"""
{whatshap} phase --distrust-genotypes --reference {ref} -o {sample}/{sample}.whatshap.vcf.gz {sample}/{sample}.vcf.gz {bam}
{bcftools} index {sample}/{sample}.whatshap.vcf.gz
{whatshap} stats --gtf {sample}/whatshap/{sample}.gtf --tsv {sample}/whatshap/{sample}.tsv  --block-list {sample}/whatshap/{sample}.blocks.txt {sample}/{sample}.whatshap.vcf.gz > {sample}/{sample}.whatshap.txt"""
	return(cmd)

def whatshap_haplotag(bam,sample,ref,samtools,whatshap):

	cmd=f"""
{whatshap} haplotag --ref {ref} -o {sample}/{sample}.happlotagged.bam {sample}/{sample}.whatshap.vcf.gz {sample}/{sample}.bam
{samtools} index {sample}/{sample}.happlotagged.bam
"""
	return(cmd)

def nanostats_qc(bam,sample,nanostats):
	cmd=f"{nanostats} --bam {bam} > {sample}/{sample}.nanostats.tsv"
	return(cmd)

def fastqc_qc(sample,fastqc):
	cmd=f"{fastqc} --memory 5000 -t 16 -o {sample} {sample}/{sample}.fastq.gz --dir {sample}/tmp/"
	return(cmd)

def multiqc_collect(sample,multiqc):
	cmd=f"{multiqc} {sample} --outdir {sample}"
	return(cmd)


def main():

	input_folder=sys.argv[1]
	sample=sys.argv[2]
	bam=f"{sample}/{sample}.bam"

	try:
		os.system(f"mkdir {sample}")
		os.system(f"mkdir {sample}/{sample}")
		os.system(f"mkdir {sample}/logs")
		os.system(f"mkdir {sample}/whatshap")
		os.system(f"mkdir {sample}/scripts")
		os.system(f"mkdir {sample}/fastq")
		os.system(f"mkdir {sample}/tmp")

	except:
		pass

	ubam_path=open(f"{sample}/bam_input.txt","w")
	ubam_path.write("\n".join(glob.glob(f"{input_folder}/*bam")))
	ubam_path.close()

	align_cmd=align(input_folder,sample,ref,samtools)

	print(align_cmd)
	align_job=Slurm("align", {"account": account, "ntasks": "20","time":"3-00:00:00" },log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	align_job_id=align_job.run(align_cmd)
	print(align_job_id)

	sniffles_cmd=sniffles_sv(bam,sample,sniffles)
	print(sniffles_cmd)
	sniffles_job=Slurm("sniffles", {"account": account, "ntasks": "16","time":"1-00:00:00","dependency":f"afterok:{align_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	sniffles_job_id=sniffles_job.run(sniffles_cmd)

	picard_cmd=picard_qc(bam,sample,ref,picard)
	print(picard_cmd)
	picard_job=Slurm("picard", {"account": account, "ntasks": "8","time":"1-00:00:00","dependency":f"afterok:{align_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	picard_job_id=picard_job.run(picard_cmd)

	pytor_cmd=cnvpytor_call(bam,sample,ref,pytor)
	print(pytor_cmd)
	pytor_job=Slurm("cnvpytor", {"account": account, "ntasks": "4","time":"1-00:00:00","dependency":f"afterok:{align_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	pytor_job_id=pytor_job.run(pytor_cmd)

	trgt_cmd=target_type(bam,sample,ref,trgt_repeats,trgt)
	print(trgt_cmd)
	trgt_job=Slurm("trgt", {"account": account, "ntasks": "1","time":"1-00:00:00","dependency":f"afterok:{align_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	trgt_job_id=trgt_job.run(trgt_cmd)

	dv_jobs=[]
	chromosomes=list(range(1,23))+["X","Y","MT"]
	for chromosome in chromosomes:
		dv_cmd=deepvariant_call_chr(bam,sample,ref,deepvariant,chromosome)
		print(dv_cmd)
		dv_job=Slurm(f"deepvariant_{chromosome}", {"account": account, "ntasks": "16","time":"1-00:00:00","dependency":f"afterok:{align_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
		dv_jobs.append(dv_job.run(dv_cmd))

	bcftools_cmd=bcftools_concat(sample,bcftools)
	print(bcftools_cmd)

	dv_jobs=":".join(map(str,dv_jobs))
	bcftools_concat_job=Slurm("concat", {"account": account, "ntasks": "1","time":"1-00:00:00","dependency":f"afterok:{dv_jobs}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	bcftools_concat_job_id=bcftools_concat_job.run(bcftools_cmd)

	whatshap_phase_cmd=whatshap_phase(bam,sample,ref,bcftools,whatshap)
	print(whatshap_phase_cmd)
	whatshap_phase_job=Slurm("whatshap_phase", {"account": account, "ntasks": "4","time":"1-00:00:00","dependency":f"afterok:{bcftools_concat_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	whatshap_phase_job_id=whatshap_phase_job.run(whatshap_phase_cmd)

	whatshap_haplotag_cmd=whatshap_haplotag(bam,sample,ref,samtools,whatshap)
	print(whatshap_haplotag_cmd)
	whatshap_haplotag_job=Slurm("whatshap_haplotag", {"account": account, "ntasks": "4","time":"1-00:00:00","dependency":f"afterok:{whatshap_phase_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	whatshap_haplotag_job_id=whatshap_haplotag_job.run(whatshap_haplotag_cmd)

	nanostats_cmd=nanostats_qc(bam,sample,nanostats)
	print(nanostats_cmd)
	nanostats_job=Slurm("nanostat", {"account": account, "ntasks": "8","time":"1-00:00:00","dependency":f"afterok:{align_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	nanostats_job_id=nanostats_job.run(nanostats_cmd)

	fastqc_qc_cmd=fastqc_qc(sample,fastqc)
	print(fastqc_qc_cmd)
	fastqc_qc_job=Slurm("fastqc", {"account": account, "ntasks": "16","time":"1-00:00:00","dependency":f"afterok:{align_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	fastqc_qc_job_id=fastqc_qc_job.run(fastqc_qc_cmd)

	multiqc_cmd=multiqc_collect(sample,multiqc)
	print(multiqc_cmd)
	multiqc_job=Slurm("multiqc", {"account": account, "ntasks": "1","time":"1-00:00:00","dependency":f"afterok:{fastqc_qc_job_id}:{whatshap_phase_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	multiqc_job_id=multiqc_job.run(multiqc_cmd)


main()
