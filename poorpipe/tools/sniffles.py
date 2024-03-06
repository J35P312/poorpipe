def sniffles_sv(bam,sample,sniffles,bcftools,wd):
	cmd=f"""
{sniffles} -m {bam} --cluster --genotype --ignore_sd --report_str -s 3 -r 500 -l 50 -t 10 -v {sample}/{sample}.snf.vcf.tmp --min_het_af 0
python {wd}/utils/clean_sniffles.py {sample}/{sample}.snf.vcf.tmp | {bcftools} reheader -s {sample}/sample_id.txt - > {sample}/{sample}.snf.vcf
"""

	return(cmd)

