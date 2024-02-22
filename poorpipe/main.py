from slurmpy import Slurm
import sys
import time
import os
import glob
import json

from tools.align import align
from tools.cnvpytor import cnvpytor_call
from tools.cnvpytor import cnvpytor_baf
from tools.sniffles import sniffles_sv

with open(sys.argv[3]) as config_file:
	config = json.load(config_file)


#extract slurm config
account=config["slurm"]["account"]
processes=config["processes"]

#setup tools
singularity_options=config["singularity"]["options"]
singularity_cmd=f"singularity exec {singularity_options}" 

samtools    =f"{singularity_cmd} { config['tools']['samtools']['singularity'] } samtools"
sniffles    =f"{singularity_cmd} { config['tools']['sniffles']['singularity'] } sniffles"
picard      =f"{singularity_cmd} { config['tools']['picard']['singularity'] } picard"
pytor       =f"{singularity_cmd} { config['tools']['cnvpytor']['singularity'] } cnvpytor"
deepvariant =f"{singularity_cmd} { config['tools']['deepvariant']['singularity'] }"
bcftools    =f"{singularity_cmd} { config['tools']['bcftools']['singularity'] } bcftools"
nanostats   =f"{singularity_cmd} { config['tools']['nanostats']['singularity'] } NanoStat"
fastqc      =f"{singularity_cmd} { config['tools']['fastqc']['singularity'] } fastqc"
multiqc     =f"{singularity_cmd} { config['tools']['multiqc']['singularity'] } multiqc"
whatshap    =f"{singularity_cmd} { config['tools']['whatshap']['singularity'] } whatshap"
bgzip       =f"{singularity_cmd} { config['tools']['bcftools']['singularity'] } bgzip"
svdb        =f"{singularity_cmd} { config['tools']['svdb']['singularity'] } svdb"
tiddit      =f"{singularity_cmd} { config['tools']['tiddit']['singularity'] } tiddit"
methylartist=f"{singularity_cmd} { config['tools']['methylartist']['singularity']} methylartist"
genmod      =f"{singularity_cmd} { config['tools']['genmod']['singularity']} genmod"
genmod_rank_model =config['tools']['genmod']["rank_model"]

trgt=config["tools"]["trgt"]["binary"]
trgt_repeats=config["tools"]["trgt"]["repeat_catalog"]

vep="{} {} vep".format( singularity_cmd,config["tools"]["vep"]["singularity"] )
filter_vep="{} {} filter_vep".format( singularity_cmd,config["tools"]["vep"]["singularity"] )
vep_options=config["tools"]["vep"]["vep_snv"]
vep_sv_options=config["tools"]["vep"]["vep_sv"]

ref=config["reference"]

def tiddit_coverage(bam,tiddit):
	cmd=f"{tiddit} --cov --bam {bam} -o {bam}\n{tiddit} --cov --bam {bam} -z 10000 -o {bam}.10kbp -w "
	return(cmd)

def methylartist_wgmeth(bam,ref,methylartist,cores):
	cmd=f"""
{methylartist} wgmeth -b {bam} -r {ref} -f {ref}.fai  --motif CG --mod m --primary_only -p {cores} -q 0 -o {bam}.m.methyl.bed
{methylartist} wgmeth -b {bam} -r {ref} -f {ref}.fai  --motif CG --mod m --primary_only -p {cores} --dss -q 0 -o {bam}.m.dss.bed
"""
	return(cmd)

def methylartist_wgmeth_phased(bam,ref,methylartist,cores):
	cmd=f"{methylartist} wgmeth -b {bam} -r {ref} -f {ref}.fai  --motif CG --mod m --primary_only -p {cores} -q 0 --phased --dss"
	return(cmd)

def target_type(bam,sample,ref,trgt_repeats,trgt,samtools):
	cmd=f"""
{samtools} view {bam} -bh -x MM -x ML >  {bam}.nometh.bam
{samtools} index {bam}.nometh.bam
{trgt} --genome {ref} --repeats {trgt_repeats} --reads {bam}.nometh.bam --output-prefix {sample}/{sample}.trgt"""
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

def deepvariant_call_chr(bam,sample,ref,deepvariant,chromosome):
	cmd=f"""
TMPDIR={sample}/tmp/
{deepvariant} run_deepvariant --output_vcf {sample}/{sample}/dv.{chromosome}.vcf --model_type ONT_R104 --num_shards 16 --reads {bam} --ref {ref} --regions {chromosome} --intermediate_results_dir {sample}/tmp/{chromosome}"""
	return(cmd)

def bcftools_concat(sample,ref,bcftools):

	bcftools_cmd=f"""
{bcftools} concat {sample}/{sample}/dv.*vcf | grep -v RefCall | {bcftools} sort -O z - --temp-dir {sample}/tmp/ | {bcftools} norm -O z -f {ref} -m- - > {sample}/{sample}.vcf.gz
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

#def fastqc_qc(sample,fastqc):
#	cmd=f"{fastqc} --memory 5000 -t 16 -o {sample} {sample}/{sample}.fastq.gz --dir {sample}/tmp/ --min_length 500"
#	return(cmd)

def multiqc_collect(sample,multiqc):
	cmd=f"{multiqc} {sample} --cl-config \"log_filesize_limit: 2000000000\" --outdir {sample} --ignore \"*tmp/\" --ignore \"*fastq/\""
	return(cmd)

def vep_annotation(sample,in_vcf,out_vcf,vep,vep_options,bcftools,bgzip):
	cmd=f"""
{vep} --offline --cache -i {in_vcf} -o {out_vcf} --force_overwrite --fork 16 --vcf {vep_options}
{bgzip} -f {out_vcf}
{bcftools} index {out_vcf}.gz

"""
	return(cmd)

def genmod_rank_snv(in_vcf,out_vcf,sample,family,filter_vep,bcftools,genmod_rank_model,genmod,bgzip):

	filter_vep_vcf=in_vcf.replace(".vcf.gz",".filter_vep.vcf")
	bcftools_quality=in_vcf.replace(".vcf.gz",".hiq.vcf.gz")
	genmod_models_vcf=in_vcf.replace(".vcf.gz",".genmod_models.vcf")

	cmd=f"""
{bcftools} view -e 'FMT/GQ < 10 | FMT/AD < 3 | INFO/AF > 0.05' {in_vcf} -o {bcftools_quality}
{genmod} annotate {bcftools_quality} --annotate_regions | {genmod} models - -f {sample}/{family}.fam  -t ped -p 10 -o {genmod_models_vcf}
{genmod} score -f {sample}/{family}.fam -t ped -c {genmod_rank_model} -r {genmod_models_vcf} -o {out_vcf}

{bgzip} -f {out_vcf}
{bcftools} index {out_vcf}.gz

"""
	return(cmd)

def svdb_query(in_vcf,out_vcf,db,svdb):
	cmd=f"{svdb} --query --query_vcf {in_vcf} --overlap 0.5 --bnd_distance 10000 --db {db} > {out_vcf}"
	return(cmd)

def submit_job(cmd,job_name,slurm_settings,sample):
	job=Slurm(job_name,slurm_settings,log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	job_id=job.run(cmd)
	return(job_id)

def get_chromosomes(fai):
	chromosomes=[]
	min_size=50000
	for line in open(fai):
		content=line.strip()

	return(chromosomes)

def finnish(sample,status,jobs):
	j=",".join(list(map(str,jobs)))

	cmd=f"""
sacct --format=jobid,jobname%50,account,partition,alloccpus,TotalCPU,elapsed,start,end,state,exitcode --jobs {j} | perl -nae 'my @headers=(jobid,jobname,account,partition,alloccpus,TotalCPU,elapsed,start,end,state,exitcode); if($. == 1) {{ print q{{#}} . join(qq{{\\t}}, @headers), qq{{\\n}} }} if ($. >= 3 && $F[0] !~ /( .batch | .bat+ )\\b/xms) {{ print join(qq{{\\t}}, @F), qq{{\\n}} }}' > {status}
"""
	return(cmd)


def main():

	input_folder=sys.argv[1]
	sample=sys.argv[2]
	family=sys.argv[4]
	bam=f"{sample}/{sample}.bam"
	jobs=[]

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

	sample_id_path=open(f"{sample}/sample_id.txt","w")
	sample_id_path.write(f"{sample}")
	sample_id_path.close()

	wd=os.path.dirname(os.path.realpath(__file__))

	readlength_script=f"{wd}/utils/seqlen.py"
	generate_ped_script=f"{wd}/utils/make_ped.py"
	align_cmd=align(input_folder,sample,family,ref,samtools,readlength_script,generate_ped_script)
	align_job_id=submit_job(align_cmd,"align",{"account": account, "ntasks":processes['align']['cores'],"time":processes['align']['time'] },sample)
	jobs.append(align_job_id)

	methylartist_wgmeth_cmd=methylartist_wgmeth(bam,ref,methylartist,16)
	methylartist_wgmeth_job_id=submit_job(methylartist_wgmeth_cmd, "methylartist_wgmeth",
	{"account": account, "ntasks": processes['methylartist_wgmeth']['cores'],"time":processes['methylartist_wgmeth']['time'],"dependency":f"afterok:{align_job_id}" },sample)
	jobs.append(methylartist_wgmeth_job_id)


	tiddit_cov_cmd=tiddit_coverage(bam,tiddit)
	tiddit_cov_job=Slurm("tiddit", {"account": account, "ntasks": "2","time":"1-00:00:00","dependency":f"afterok:{align_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	tiddit_cov_job_id=tiddit_cov_job.run(tiddit_cov_cmd)
	jobs.append(tiddit_cov_job_id)

	sniffles_cmd=sniffles_sv(bam,sample,sniffles,bcftools)
	print(sniffles_cmd)
	sniffles_job=Slurm("sniffles", {"account": account, "ntasks": "16","time":"1-00:00:00","dependency":f"afterok:{align_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	sniffles_job_id=sniffles_job.run(sniffles_cmd)
	jobs.append(sniffles_job_id)

	picard_cmd=picard_qc(bam,sample,ref,picard)
	print(picard_cmd)
	picard_job=Slurm("picard", {"account": account, "ntasks": "8","time":"1-00:00:00","dependency":f"afterok:{align_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	picard_job_id=picard_job.run(picard_cmd)
	jobs.append(picard_job_id)

	picard_gc_cmd=picard_gc(bam,sample,ref,picard,samtools)
	picard_gc_job=Slurm("picard_gc", {"account": account, "ntasks": "8","time":"1-00:00:00","dependency":f"afterok:{align_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	picard_gc_job_id=picard_gc_job.run(picard_gc_cmd)
	jobs.append(picard_gc_job_id)

	pytor_cmd=cnvpytor_call(bam,sample,ref,pytor)
	print(pytor_cmd)
	pytor_job=Slurm("cnvpytor", {"account": account, "ntasks": "4","time":"1-00:00:00","dependency":f"afterok:{align_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	pytor_job_id=pytor_job.run(pytor_cmd)
	jobs.append(pytor_job_id)


	trgt_cmd=target_type(bam,sample,ref,trgt_repeats,trgt,samtools)
	print(trgt_cmd)
	trgt_job=Slurm("trgt", {"account": account, "ntasks": "1","time":"1-00:00:00","dependency":f"afterok:{align_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	trgt_job_id=trgt_job.run(trgt_cmd)
	jobs.append(trgt_job_id)


	dv_jobs=[]
	chromosomes=list(range(1,23))+["X","Y","MT"]
	for chromosome in chromosomes:
		dv_cmd=deepvariant_call_chr(bam,sample,ref,deepvariant,chromosome)
		print(dv_cmd)
		dv_job=Slurm(f"deepvariant_{chromosome}", {"account": account, "ntasks": "20","time":"1-00:00:00","dependency":f"afterok:{align_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
		dv_jobs.append(dv_job.run(dv_cmd))
	jobs+=dv_jobs

	bcftools_cmd=bcftools_concat(sample,ref,bcftools)
	print(bcftools_cmd)

	dv_jobs=":".join(map(str,dv_jobs))
	bcftools_concat_job=Slurm("concat", {"account": account, "ntasks": "1","time":"1-00:00:00","dependency":f"afterok:{dv_jobs}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	bcftools_concat_job_id=bcftools_concat_job.run(bcftools_cmd)
	jobs.append(bcftools_concat_job_id)


	pytor_baf_cmd=cnvpytor_baf(bam,sample,ref,pytor)
	pytor_baf_job=Slurm("cnvpytor_baf", {"account": account, "ntasks": "1","time":"1-00:00:00","dependency":f"afterok:{bcftools_concat_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	pytor_baf_job_id=pytor_baf_job.run(pytor_baf_cmd)
	jobs.append(pytor_baf_job_id)


	svdb_merge_cmd=svdb_merge(f"{sample}/{sample}.pytor.vcf",f"{sample}/{sample}.snf.vcf",svdb,sample)
	svdb_merge_job=Slurm("svdb_merge", {"account": account, "ntasks": "1","time":"1-00:00:00","dependency":f"afterok:{pytor_job_id}:{sniffles_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	svdb_merge_job_id=svdb_merge_job.run(svdb_merge_cmd)
	jobs.append(svdb_merge_job_id)


	vep_sv_cmd=vep_annotation(sample,f"{sample}/{sample}.sv.vcf",f"{sample}/{sample}.sv.vep.vcf",vep,vep_sv_options,bcftools,bgzip)
	vep_sv_job=Slurm("vep_sv_annotation", {"account": account, "ntasks": "16","time":"1-00:00:00","dependency":f"afterok:{svdb_merge_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	vep_sv_job_id=vep_sv_job.run(vep_sv_cmd)
	jobs.append(vep_sv_job_id)

	if not "" == config['tools']['svdb']['database']:
		svdb_query_cmd=svdb_query(f"{sample}/{sample}.sv.vep.vcf.gz",f"{sample}/{sample}.sv.vep.svdb.vcf",config['tools']['svdb']['database'],svdb)
		svdb_query_job=Slurm("svdb_query", {"account": account, "ntasks": "2","time":"1-00:00:00","dependency":f"afterok:{vep_sv_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
		svdb_query_job_id=svdb_query_job.run(svdb_query_cmd)
		jobs.append(svdb_query_job_id)

	whatshap_phase_cmd=whatshap_phase(bam,sample,ref,bcftools,whatshap)
	print(whatshap_phase_cmd)
	whatshap_phase_job=Slurm("whatshap_phase", {"account": account, "ntasks": "4","time":"1-00:00:00","dependency":f"afterok:{bcftools_concat_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	whatshap_phase_job_id=whatshap_phase_job.run(whatshap_phase_cmd)
	jobs.append(whatshap_phase_job_id)

	whatshap_haplotag_cmd=whatshap_haplotag(bam,sample,ref,samtools,whatshap)
	print(whatshap_haplotag_cmd)
	whatshap_haplotag_job=Slurm("whatshap_haplotag", {"account": account, "ntasks": "4","time":"1-00:00:00","dependency":f"afterok:{whatshap_phase_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	whatshap_haplotag_job_id=whatshap_haplotag_job.run(whatshap_haplotag_cmd)
	jobs.append(whatshap_haplotag_job_id)

	methylartist_wgmeth_phased_cmd=methylartist_wgmeth_phased(f"{sample}/{sample}.happlotagged.bam",ref,methylartist,16)
	methylartist_wgmeth_phased_job_id=submit_job(methylartist_wgmeth_phased_cmd, "methylartist_wgmeth",
	{"account": account, "ntasks": processes['methylartist_wgmeth']['cores'],"time":processes['methylartist_wgmeth']['time'],"dependency":f"afterok:{whatshap_haplotag_job_id}" },sample)
	jobs.append(methylartist_wgmeth_phased_job_id)

	vep_cmd=vep_annotation(sample,f"{sample}/{sample}.whatshap.vcf.gz",f"{sample}/{sample}.vcf",vep,vep_options,bcftools,bgzip)
	vep_job=Slurm("vep_annotation", {"account": account, "ntasks": "16","time":"1-00:00:00","dependency":f"afterok:{whatshap_phase_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	vep_job_id=vep_job.run(vep_cmd)
	jobs.append(vep_job_id)

	genmod_cmd=genmod_rank_snv(f"{sample}/{sample}.vcf.gz",f"{sample}/{sample}.scored.vcf",sample,family,filter_vep,bcftools,genmod_rank_model,genmod,bgzip)
	genmod_job=Slurm("genmod_annotation", {"account": account, "ntasks": "2","time":"1-00:00:00","dependency":f"afterok:{vep_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	genmod_job_id=genmod_job.run(genmod_cmd)
	jobs.append(genmod_job_id)


	nanostats_cmd=nanostats_qc(bam,sample,nanostats)
	print(nanostats_cmd)
	nanostats_job=Slurm("nanostat", {"account": account, "ntasks": "8","time":"1-00:00:00","dependency":f"afterok:{align_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	nanostats_job_id=nanostats_job.run(nanostats_cmd)
	jobs.append(nanostats_job_id)

	multiqc_cmd=multiqc_collect(sample,multiqc)
	print(multiqc_cmd)
	multiqc_job=Slurm("multiqc", {"account": account, "ntasks": "1","time":"1-00:00:00","dependency":f"afterok:{whatshap_phase_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	multiqc_job_id=multiqc_job.run(multiqc_cmd)
	jobs.append(multiqc_job_id)

	complete_cmd=finnish(sample,f"{sample}/complete",jobs)
	complete_job=Slurm(f"complete-{sample}", {"account": account, "ntasks": "1","time":"1-00:00:00","dependency":f"afterok:{multiqc_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	complete_job_id=complete_job.run(complete_cmd)

	fail_cmd=finnish(sample,f"{sample}/fail",jobs)
	fail_job=Slurm(f"fail-{sample}", {"account": account, "ntasks": "1","time":"1-00:00:00","dependency":f"afternotok:{multiqc_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	fail_job_id=fail_job.run(fail_cmd)

	time.sleep(10)
	os.system(finnish(sample,f"{sample}/submitted",jobs))

main()
