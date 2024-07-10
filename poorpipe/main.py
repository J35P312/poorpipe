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
cytosure    =f"{singularity_cmd} { config['tools']['vcf2cytosure']['singularity'] } vcf2cytosure"
clair3 =f"{singularity_cmd} { config['tools']['clair3']['singularity'] } run_clair3.sh "
hificnv = f"{singularity_cmd} { config['tools']['hificnv']['singularity'] } hificnv"
bcftools    =f"{singularity_cmd} { config['tools']['bcftools']['singularity'] } bcftools"
nanostats   =f"{singularity_cmd} { config['tools']['nanostats']['singularity'] } NanoStat"
fastqc      =f"{singularity_cmd} { config['tools']['fastqc']['singularity'] } fastqc"
multiqc     =f"{singularity_cmd} { config['tools']['multiqc']['singularity'] } multiqc"
whatshap    =f"{singularity_cmd} { config['tools']['whatshap']['singularity'] } whatshap"
bgzip       =f"{singularity_cmd} { config['tools']['bcftools']['singularity'] } bgzip"
svdb        =f"{singularity_cmd} { config['tools']['svdb']['singularity'] } svdb"
tiddit      =f"{singularity_cmd} { config['tools']['tiddit']['singularity'] } tiddit"
chopper      =f"{singularity_cmd} { config['tools']['chopper']['singularity'] } chopper"
methylartist=f"{singularity_cmd} { config['tools']['methylartist']['singularity']} methylartist"

paraphase=f"{singularity_cmd} { config['tools']['paraphase']['singularity']} paraphase"
paraphase_genome=config['tools']['paraphase']['genome']

cn_male =   config['tools']['hificnv']['cn_male']
cn_female = config['tools']['hificnv']['cn_female']


minimap2    =f"{singularity_cmd} { config['tools']['minimap2']['singularity']} minimap2"
vcfsort     =f"{singularity_cmd} { config['tools']['vcftools']['singularity']} vcf-sort"
genmod      =f"{singularity_cmd} { config['tools']['genmod']['singularity']} genmod"
genmod_rank_model =config['tools']['genmod']["rank_model"]
genmod_sv_rank_model =config['tools']['genmod']["sv_rank_model"]

strdust=f"{singularity_cmd} { config['tools']['strdust']['singularity'] } STRdust"
strdust_repeats=config["tools"]["strdust"]["repeat_catalog"]

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

def paraphase_type(sample,ref,paraphase_genome,paraphase):
	cmd=f"{paraphase} -b {sample}/{sample}.bam -r {ref} -o {sample}/paraphase --genome {paraphase_genome}"
	return(cmd)

def run_hificnv(sample,cn_male,cn_female,ref,hificnv):
	cmd=f"""
if grep -q male {sample}/{sample}.gender.txt; then
	{hificnv} --bam {sample}/{sample}.bam --ref {ref} --output-prefix {sample}/{sample}.hificnv --threads 8 --maf {sample}/{sample}.pythor_snp.vcf.gz  --expected-cn {cn_male}
else
	{hificnv} --bam {sample}/{sample}.bam --ref {ref} --output-prefix {sample}/{sample}.hificnv --threads 8 --maf {sample}/{sample}.pythor_snp.vcf.gz  --expected-cn {cn_female}
fi


"""
	return(cmd)

def cytosure_cgh(sample,coverage_bed,sv_vcf,cytosure):

	pytor_vcf=sv_vcf.replace(".vcf.gz","pytor.vcf")

	cmd=f"""
zgrep -E 'pytor|#' {sv_vcf} > {pytor_vcf}
sex==$(cat {sample}/{sample}.gender.txt)
{cytosure} --size 5000 --frequency 0.1 --coverage {coverage_bed} --vcf {pytor_vcf} --genome 37 --bins 40 --sex $sex
"""
	return(cmd)

def target_type(bam,sample,ref,strdust_repeats,strdust):

	cmd=f"{strdust} --region-file {strdust_repeats} {ref} {bam} --sample {sample} > {sample}/{sample}.strdust.vcf"

	return(cmd)

def picard_qc(bam,sample,ref,picard):
	cmd=f"{picard}  -Xmx20G CollectWgsMetrics --INPUT {bam} --OUTPUT {sample}/{sample}.wgsmetrics.txt --REFERENCE_SEQUENCE {ref} --COUNT_UNPAIRED true --MINIMUM_BASE_QUALITY 1 --VALIDATION_STRINGENCY SILENT"
	return(cmd)

def svdb_merge(pytor,sniffles,svdb,sample):
	cmd=f"{svdb} --merge --vcf {pytor}:cnvpytor {sniffles}:sniffles --priority sniffles,cnvpytor --no_var --no_intra --bnd_distance 5000 --overlap 0.5 > {sample}/{sample}.sv.vcf"
	return(cmd)

def picard_gc(bam,sample,ref,picard,samtools):
	cmd=f"""
{samtools} view -bh -F 256 {bam} > {sample}/{sample}.filt.bam
{samtools} index {sample}/{sample}.filt.bam
{picard} -Xmx20G CollectGcBiasMetrics --CHART {sample}/{sample}.gc.pdf --INPUT {sample}/{sample}.filt.bam --SUMMARY_OUTPUT {sample}/{sample}.gc_summary.txt --OUTPUT {sample}/{sample}.gc.txt --REFERENCE_SEQUENCE {ref} --VALIDATION_STRINGENCY SILENT
"""
	return(cmd)

def bcftools_concat(sample,ref,bcftools):

	bcftools_cmd=f"""
{bcftools} concat {sample}/{sample}_clair3/merge_output.vcf.gz | {bcftools} sort -O z - --temp-dir {sample}/tmp/ | {bcftools} norm -O z -f {ref} -m- - > {sample}/{sample}.vcf.gz
{bcftools} index {sample}/{sample}.vcf.gz
{bcftools} view  {sample}/{sample}.vcf.gz -i "QUAL>20 & FORMAT/DP > 20 & FORMAT/GQ > 10" > {sample}/{sample}.pytor.vcf.gz

{bcftools} view {sample}/{sample}.vcf.gz > {sample}/{sample}.pythor_snp.vcf.gz
{bcftools} stats {sample}/{sample}.vcf.gz > {sample}/{sample}.txt
	
"""
	return(bcftools_cmd)

def run_clair3(bam,sample,ref,model,platform,clair3):
	cmd=f"""
{clair3} --output={sample}/{sample}_clair3 --sample_name={sample} --bam_fn={bam} --ref_fn {ref} --threads 18 -m {model} --platform {platform}
	"""
	return(cmd)

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
	cmd=f"{nanostats} --fastq {sample}/{sample}.fastq.gz > {sample}/{sample}.nanostats.tsv"
	return(cmd)

def multiqc_collect(sample,multiqc):
	cmd=f"{multiqc} {sample} --cl-config \"log_filesize_limit: 2000000000\" --outdir {sample} --ignore \"*tmp/\" --ignore \"*fastq/\""
	return(cmd)

def vep_annotation(sample,in_vcf,out_vcf,vep,vep_options,bcftools,bgzip):
	cmd=f"""
{vep} --offline --cache -i {in_vcf} -o {out_vcf} --force_overwrite --fork 16 --vcf {vep_options} --format vcf
{bgzip} -f {out_vcf}
{bcftools} index {out_vcf}.gz

"""
	return(cmd)


def genmod_rank_snv(in_vcf,out_vcf,sample,family,filter_vep,bcftools,genmod_rank_model,genmod,bgzip):

	filter_vep_vcf=in_vcf.replace(".vcf.gz",".filter_vep.vcf")
	bcftools_quality=in_vcf.replace(".vcf.gz",".hiq.vcf.gz")
	genmod_models_vcf=in_vcf.replace(".vcf.gz",".genmod_models.vcf")

	cmd=f"""
{bcftools} view -e 'FMT/GQ < 5 | FMT/DP < 5 | FORMAT/AF < 0.05' {in_vcf} | {bcftools} +split-vep -c AF:Float -e 'AF>0.2' -O z - > {bcftools_quality}
{genmod} annotate {bcftools_quality} --annotate_regions | {genmod} models - -f {sample}/{family}.fam  -t ped -p 10 -o {genmod_models_vcf}
{genmod} score -f {sample}/{family}.fam -t ped -c {genmod_rank_model} -r {genmod_models_vcf} -o {out_vcf}

{bgzip} -f {out_vcf}
{bcftools} index {out_vcf}.gz

"""
	return(cmd)

def genmod_rank_sv(in_vcf,out_vcf,sample,family,filter_vep,bcftools,genmod_rank_model,genmod,bgzip,vcfsort):

	bcftools_quality=in_vcf.replace(".sv.vep.svdb.vcf",".sv.hiq.vcf.gz")
	genmod_models_vcf=in_vcf.replace(".sv.vep.svdb.vcf",".sv.genmod_models.vcf")

	cmd=f"""
{bcftools} view -e 'INFO/FRQ > 0.1' {in_vcf} -o {bcftools_quality}
{genmod} annotate {bcftools_quality} --annotate_regions | {genmod} models - -f {sample}/{family}.fam  -t ped -p 10 -o {genmod_models_vcf}
{genmod} score -f {sample}/{family}.fam -t ped -c {genmod_rank_model} -r {genmod_models_vcf} -o {out_vcf}

{vcfsort} {out_vcf} > {out_vcf}.tmp
mv {out_vcf}.tmp {out_vcf}
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

def finnish(sample,to_remove,status,jobs):
	j=",".join(list(map(str,jobs)))

	rm="\nrm -rf ".join(to_remove)
	if to_remove:
		rm+= "rm -rf "

	cmd=f"""
sacct --format=jobid,jobname%50,account,partition,alloccpus,TotalCPU,elapsed,start,end,state,exitcode --jobs {j} | perl -nae 'my @headers=(jobid,jobname,account,partition,alloccpus,TotalCPU,elapsed,start,end,state,exitcode); if($. == 1) {{ print q{{#}} . join(qq{{\\t}}, @headers), qq{{\\n}} }} if ($. >= 3 && $F[0] !~ /( .batch | .bat+ )\\b/xms) {{ print join(qq{{\\t}}, @F), qq{{\\n}} }}' > {status}
{rm}
"""
	return(cmd)


def main():

	input_folder=sys.argv[1]
	sample=sys.argv[2]
	family=sys.argv[4]
	if len(sys.argv) > 5:
		method=sys.argv[5]
	else:
		method="ont"
	if method != "ont" and method != "hifi":
		print("ont and hifi supported")
		quit()
	

	bam=f"{sample}/{sample}.bam"
	jobs=[]

	try:
		os.system(f"mkdir {sample}")
		os.system(f"mkdir {sample}/{sample}")
		os.system(f"mkdir {sample}/logs")
		try:
			os.system(f"rm {sample}/logs/*")
		except:
			pass


		os.system(f"mkdir {sample}/whatshap")
		os.system(f"mkdir {sample}/scripts")

		try:
			os.system(f"rm {sample}/scripts/*")
		except:
			pass

		os.system(f"mkdir {sample}/fastq")
		os.system(f"mkdir {sample}/tmp")
		os.system(f"mkdir {sample}/pytor")
		try:
			os.system(f"rm {sample}/fail")
		except:
			pass
		try:
			os.system(f"rm {sample}/complete")
		except:
			pass


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
	align_cmd=align(input_folder,sample,family,ref,samtools,readlength_script,generate_ped_script,minimap2,chopper,method)
	align_job_id=submit_job(align_cmd,"align",{"account": account, "ntasks":processes['align']['cores'],"time":processes['align']['time'] },sample)
	jobs.append(align_job_id)

	paraphase_type_cmd=paraphase_type(sample,ref,paraphase_genome,paraphase)
	paraphase_type_job_id=submit_job(paraphase_type_cmd, "paraphase_type",
	{"account": account, "ntasks": "1","time":"1-00:00:00","dependency":f"afterok:{align_job_id}" },sample)
	jobs.append(paraphase_type_job_id)

	methylartist_wgmeth_cmd=methylartist_wgmeth(bam,ref,methylartist,16)
	methylartist_wgmeth_job_id=submit_job(methylartist_wgmeth_cmd, "methylartist_wgmeth",
	{"account": account, "ntasks": processes['methylartist_wgmeth']['cores'],"time":processes['methylartist_wgmeth']['time'],"dependency":f"afterok:{align_job_id}" },sample)
	jobs.append(methylartist_wgmeth_job_id)


	tiddit_cov_cmd=tiddit_coverage(bam,tiddit)
	tiddit_cov_job=Slurm("tiddit", {"account": account, "ntasks": "2","time":"1-00:00:00","dependency":f"afterok:{align_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	tiddit_cov_job_id=tiddit_cov_job.run(tiddit_cov_cmd)
	jobs.append(tiddit_cov_job_id)

	sniffles_cmd=sniffles_sv(bam,sample,sniffles,bcftools,wd)
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

	if method == "ont":
		clair3_platform="ont"
		clair3_model=config['tools']['clair3']['ont_model']
	else:
		clair3_platform="hifi"
		clair3_model="/opt/models/hifi_revio"

	clair3_cmd=run_clair3(bam,sample,ref,clair3_model,clair3_platform,clair3)
	clair3_job=Slurm(f"clair3", {"account": account, "ntasks": "18","time":"2-00:00:00","dependency":f"afterok:{align_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	clair3_job_id=clair3_job.run(clair3_cmd)
	jobs.append(clair3_job_id)

	bcftools_cmd=bcftools_concat(sample,ref,bcftools)
	print(bcftools_cmd)

	bcftools_concat_job=Slurm("concat", {"account": account, "ntasks": "1","time":"1-00:00:00","dependency":f"afterok:{clair3_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	bcftools_concat_job_id=bcftools_concat_job.run(bcftools_cmd)
	jobs.append(bcftools_concat_job_id)

	pytor_baf_cmd=cnvpytor_baf(bam,sample,ref,pytor)
	pytor_baf_job=Slurm("cnvpytor_baf", {"account": account, "ntasks": "1","time":"1-00:00:00","dependency":f"afterok:{bcftools_concat_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	pytor_baf_job_id=pytor_baf_job.run(pytor_baf_cmd)
	jobs.append(pytor_baf_job_id)

	hificnv_cmd=run_hificnv(sample,cn_male,cn_female,ref,hificnv)
	hificnv_job=Slurm("hificnv", {"account": account, "ntasks": "1","time":"1-00:00:00","dependency":f"afterok:{bcftools_concat_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	hificnv_job_id=hificnv_job.run(hificnv_cmd)
	jobs.append(hificnv_job_id)

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

		genmod_sv_cmd=genmod_rank_sv(f"{sample}/{sample}.sv.vep.svdb.vcf",f"{sample}/{sample}.sv.scored.vcf",sample,family,filter_vep,bcftools,genmod_sv_rank_model,genmod,bgzip,vcfsort)
		genmod_sv_job=Slurm("genmod_sv_annotation", {"account": account, "ntasks": "2","time":"1-00:00:00","dependency":f"afterok:{svdb_query_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
		genmod_sv_job_id=genmod_sv_job.run(genmod_sv_cmd)
		jobs.append(genmod_sv_job_id)

		cytosure_cmd=cytosure_cgh(sample,f"{bam}.bed",f"{sample}/{sample}.sv.scored.vcf.gz",cytosure)
		cytosure_job=Slurm("cytosure_cgh", {"account": account, "ntasks": "1","time":"1-00:00:00","dependency":f"afterok:{genmod_sv_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
		cytosure_job_id=cytosure_job.run(cytosure_cmd)
		jobs.append(cytosure_job_id)


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

	strdust_cmd=target_type(f"{sample}/{sample}.happlotagged.bam",sample,ref,strdust_repeats,strdust)
	print(strdust_cmd)
	strdust_job=Slurm("strdust", {"account": account, "ntasks": "1","time":"1-00:00:00","dependency":f"afterok:{whatshap_haplotag_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	strdust_job_id=strdust_job.run(strdust_cmd)
	jobs.append(strdust_job_id)

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
	multiqc_job=Slurm("multiqc", {"account": account, "ntasks": "1","time":"1-00:00:00","dependency":f"afterok:{whatshap_phase_job_id}:{nanostats_job_id}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	multiqc_job_id=multiqc_job.run(multiqc_cmd)
	jobs.append(multiqc_job_id)

	to_remove=[f"{sample}/tmp",f"{sample}/fastq",f"{sample}/{sample}.fastq.gz",f"{sample}/{sample}.bam.nometh.bam"]

	all_jobs=":".join(list(map(str,jobs)))
	complete_cmd=finnish(sample,to_remove,f"{sample}/complete",jobs)
	complete_job=Slurm(f"complete-{sample}", {"account": account, "ntasks": "1","time":"00:10:00","dependency":f"afterok:{all_jobs}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	complete_job_id=complete_job.run(complete_cmd)

	fail_cmd=finnish(sample,to_remove,f"{sample}/fail",jobs)
	fail_job=Slurm(f"fail-{sample}", {"account": account, "ntasks": "1","time":"00:10:00","dependency":f"afternotok:{all_jobs}"},log_dir=f"{sample}/logs",scripts_dir=f"{sample}/scripts")
	fail_job_id=fail_job.run(fail_cmd)

	time.sleep(10)
	os.system(finnish(sample,[],f"{sample}/submitted",jobs))

main()
