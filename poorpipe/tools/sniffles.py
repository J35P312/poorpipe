def sniffles_sv(bam,sample,sniffles):
	#cmd=f"""{sniffles} -i {bam} --snf {sample}/{sample}.snf --vcf {sample}/{sample}.snf.vcf --minsupport 3 --max-del-seq-len 100 --allow-overwrite  --qc-stdev False --long-ins-length 50000 --long-del-length 999999999 --long-dup-length 999999999 --bnd-min-split-length 500 --min-alignment-length 500  --sample-id {sample}"""
	cmd=f"""{sniffles} -m {bam} --cluster --genotype --ignore_sd --report_str -s 3 -r 500 -l 50 -t 10 -v {sample}/{sample}.snf --min_het_af 0"""

	return(cmd)
