from slurmpy import Slurm
import sys
import os
import glob
import json

with open(sys.argv[3]) as config_file:
	config = json.load(config_file)


#extract slurm config
account=config["slurm"]["account"]

#setup tools
singularity_options=config["singularity"]["options"]
singularity_cmd=f"singularity exec {singularity_options}" 

samtools=f"{singularity_cmd} singularity/samtools_1.19--h50ea8bc_0.sif samtools"
sniffles=f"{singularity_cmd} singularity/sniffles_2.2--pyhdfd78af_0.sif sniffles"
picard=f"{singularity_cmd} singularity/picard_3.1.1--hdfd78af_0.sif picard"
pytor=f"{singularity_cmd} singularity/cnvpytor_1.3.1--pyhdfd78af_1.sif cnvpytor"
deepvariant=f"{singularity_cmd} singularity/deepvariant_latest.sif"
bcftools=f"{singularity_cmd} singularity/bcftools_1.19--h8b25389_0.sif bcftools"
nanostats=f"{singularity_cmd} singularity/nanostat_1.6.0--pyhdfd78af_0.sif NanoStat"
fastqc=f"{singularity_cmd} singularity/fastqc_0.12.1--hdfd78af_0.sif fastqc"
multiqc=f"{singularity_cmd} singularity/multiqc_1.18--pyhdfd78af_0.sif multiqc"
whatshap=f"{singularity_cmd} singularity/whatshap_2.1--py39h1f90b4d_0.sif whatshap"
bgzip=f"{singularity_cmd} singularity/bcftools_1.19--h8b25389_0.sif bgzip"
svdb=f"{singularity_cmd} singularity/svdb_2.8.2--py38h24c8ff8_1.sif svdb"
tiddit=f"{singularity_cmd} singularity/tiddit_3.6.1--py38h24c8ff8_0.sif tiddit"

trgt="{} {} trgt".format(singularity_cmd,config["tools"]["trgt"]["singularity"] )
trgt_repeats=config["tools"]["trgt"]["repeat_catalog"]

vep="{} {} vep".format( singularity_cmd,config["tools"]["vep"]["singularity"] )
vep_options=config["tools"]["vep"]["vep_snv"]
vep_sv_options=config["tools"]["vep"]["vep_sv"]

ref=config["reference"]

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

parallel --tmpdir {sample}/tmp -j 20 < {prefix}/parallel_fastq.txt
cat {sample}/fastq/* > {prefix}/{prefix}.fastq.gz
minimap2 -R "@RG\\tID:{prefix}\\tSM:{prefix}" -a -t 16 --MD -x map-ont {ref} {prefix}/{prefix}.fastq.gz | {samtools} sort -T {sample}/tmp/{sample}.samtools -m 5G -@8 - > {prefix}/{prefix}.bam
{samtools} index -@16  {prefix}/{prefix}.bam
{samtools} stats {prefix}/{prefix}.bam > {prefix}/{prefix}.bam.stats.txt
"""
	return(cmd)

def sniffles_sv(bam,sample,sniffles):
	cmd=f"""{sniffles} -i {bam} --snf {sample}/{sample}.snf --vcf {sample}/{sample}.snf.vcf --minsupport 3 --max-del-seq-len 100 --allow-overwrite  --qc-stdev False --long-ins-length 50000 --long-del-length 999999999 --long-dup-length 999999999 --bnd-min-split-length 500 --min-alignment-length 500  --sample-id {sample}"""

	return(cmd)

def tiddit_coverage(bam,tiddit):
	cmd=f"{tiddit} --cov --bam {bam} -o {bam}"
	return(cmd)

def target_type(bam,sample,ref,trgt_repeats,trgt):
	cmd=f"{trgt} --genome {ref} --repeats {trgt_repeats} --reads {bam} --output-prefix {sample}/{sample}.trgt"
	return(cmd)

def picard_qc(bam,sample,ref,picard):
	cmd=f"{picard}  -Xmx20G CollectWgsMetrics --INPUT {bam} --OUTPUT {sample}/{sample}.wgsmetrics.txt --REFERENCE_SEQUENCE {ref} --COUNT_UNPAIRED true --MINIMUM_BASE_QUALITY 1 --VALIDATION_STRINGENCY SILENT"
	return(cmd)

def svdb_merge(pytor,sniffles,svdb,sample):
	cmd=f"{svdb} --merge --vcf {pytor} {sniffles} --no_var --no_intra --bnd_distance 5000 --overlap 0.5 > {sample}/{sample}.sv.vcf"
	return(cmd)

def picard_gc(bam,sample,ref,picard,samtools):
	cmd=f"""
{samtools} view -bh -F 256 {bam} > {sample}/{sample}.filt.bam
{samtools} index {sample}/{sample}.filt.bam
{picard} -Xmx20G CollectGcBiasMetrics --CHART {sample}/{sample}.gc.pdf --INPUT {sample}/{sample}.filt.bam --SUMMARY_OUTPUT {sample}/{sample}.gc_summary.txt --OUTPUT {sample}/{sample}.gc.txt --REFERENCE_SEQUENCE {ref} --VALIDATION_STRINGENCY SILENT
"""
	return(cmd)

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


def deepvariant_call_chr(bam,sample,ref,deepvariant,chromosome):
	cmd=f"""
TMPDIR={sample}/tmp/
{deepvariant} run_deepvariant --output_vcf {sample}/{sample}/dv.{chromosome}.vcf --model_type ONT_R104 --num_shards 16 --reads {bam} --ref {ref} --regions {chromosome} --intermediate_results_dir {sample}/tmp/{chromosome}"""
	return(cmd)

def bcftools_concat(sample,bcftools):

	bcftools_cmd=f"""
{bcftools} concat {sample}/{sample}/dv.*vcf | grep -v RefCall | {bcftools} sort -O z - --temp-dir {sample}/tmp/ > {sample}/{sample}.vcf.gz
{bcftools} index {sample}/{sample}.vcf.gz
{bcftools} view  {sample}/{sample}.vcf.gz -i "QUAL>20 & FORMAT/DP > 20 & FORMAT/GQ > 10" > {sample}/{sample}.pytor.vcf.gz

{bcftools} view {sample}/{sample}.vcf.gz > {sample}/{sample}.pythor_snp.vcf.gz
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
	cmd=f"{multiqc} {sample} --cl-config \"log_filesize_limit: 2000000000\" --outdir {sample}"
	return(cmd)

def vep_annotation(sample,in_vcf,out_vcf,vep,vep_options,bcftools,bgzip):
	cmd=f"""
{vep} --offline --cache -i {in_vcf} -o {out_vcf} --force_overwrite --fork 16 --vcf {vep_options}
{bgzip} -f {out_vcf}
{bcftools} index {out_vcf}.gz

"""
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
		os.system(f"mkdir {sample}/pytor")

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

	tiddit_cov_cmd=tiddit_coverage(bam,tiddit)
	tiddit_cov_job=Slurm("tiddit", {"account": account, "ntasks": "2","time":"1-00:00:00","dependency":f"afterok:{align_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	tiddit_cov_job_id=tiddit_cov_job.run(tiddit_cov_cmd)

	sniffles_cmd=sniffles_sv(bam,sample,sniffles)
	print(sniffles_cmd)
	sniffles_job=Slurm("sniffles", {"account": account, "ntasks": "16","time":"1-00:00:00","dependency":f"afterok:{align_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	sniffles_job_id=sniffles_job.run(sniffles_cmd)

	picard_cmd=picard_qc(bam,sample,ref,picard)
	print(picard_cmd)
	picard_job=Slurm("picard", {"account": account, "ntasks": "8","time":"1-00:00:00","dependency":f"afterok:{align_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	picard_job_id=picard_job.run(picard_cmd)

	picard_gc_cmd=picard_gc(bam,sample,ref,picard,samtools)
	picard_gc_job=Slurm("picard_gc", {"account": account, "ntasks": "8","time":"1-00:00:00","dependency":f"afterok:{align_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	picard_gc_job_id=picard_gc_job.run(picard_gc_cmd)

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

	pytor_baf_cmd=cnvpytor_baf(bam,sample,ref,pytor)
	pytor_baf_job=Slurm("cnvpytor_baf", {"account": account, "ntasks": "1","time":"1-00:00:00","dependency":f"afterok:{bcftools_concat_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	pytor_baf_job_id=pytor_baf_job.run(pytor_baf_cmd)

	svdb_merge_cmd=svdb_merge(f"{sample}/{sample}.pytor.vcf",f"{sample}/{sample}.snf.vcf",svdb,sample)
	svdb_merge_job=Slurm("svdb_merge", {"account": account, "ntasks": "1","time":"1-00:00:00","dependency":f"afterok:{pytor_job_id}:{sniffles_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	svdb_merge_job_id=svdb_merge_job.run(svdb_merge_cmd)

	vep_sv_cmd=vep_annotation(sample,f"{sample}/{sample}.sv.vcf",f"{sample}/{sample}.sv.vep.vcf",vep,vep_sv_options,bcftools,bgzip)
	vep_sv_job=Slurm("vep_sv_annotation", {"account": account, "ntasks": "16","time":"1-00:00:00","dependency":f"afterok:{svdb_merge_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	vep_sv_job_id=vep_sv_job.run(vep_sv_cmd)

	whatshap_phase_cmd=whatshap_phase(bam,sample,ref,bcftools,whatshap)
	print(whatshap_phase_cmd)
	whatshap_phase_job=Slurm("whatshap_phase", {"account": account, "ntasks": "4","time":"1-00:00:00","dependency":f"afterok:{bcftools_concat_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	whatshap_phase_job_id=whatshap_phase_job.run(whatshap_phase_cmd)

	whatshap_haplotag_cmd=whatshap_haplotag(bam,sample,ref,samtools,whatshap)
	print(whatshap_haplotag_cmd)
	whatshap_haplotag_job=Slurm("whatshap_haplotag", {"account": account, "ntasks": "4","time":"1-00:00:00","dependency":f"afterok:{whatshap_phase_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	whatshap_haplotag_job_id=whatshap_haplotag_job.run(whatshap_haplotag_cmd)

	vep_cmd=vep_annotation(sample,f"{sample}/{sample}.whatshap.vcf.gz",f"{sample}/{sample}.vcf",vep,vep_options,bcftools,bgzip)
	vep_job=Slurm("vep_annotation", {"account": account, "ntasks": "16","time":"1-00:00:00","dependency":f"afterok:{whatshap_phase_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	vep_job_id=vep_job.run(vep_cmd)

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
