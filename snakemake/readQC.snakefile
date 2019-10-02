# Pipeline for collecting insert size, rRNA count, read genomic distribution and gene body coverage
# Command for dry run: snakemake --snakefile readQC.snakefile -n -p -j 20 --forceall
# snakemake --snakefile readQC.snakefile --dag | dot -Tpdf > dag.pdf
# snakemake --snakefile readQC.snakefile --rulegraph | dot -Tpdf > ruleGraph.pdf
# Command to execute: snakemake --snakefile readQC.snakefile -j 20
# Remember -j will replace the specific threads defined in the rules!

os.system("mkdir -p rseqcOutput")

sample = ["SRR6664591", "SRR6664592", "SRR6664615", "SRR6664616"]

rule all:
    input:
	    expand("rseqcOutput/{sample}.rRNA.txt", sample = sample),
	    expand("rseqcOutput/{sample}.readDist.txt", sample = sample),
	    expand("rseqcOutput/{sample}.size.pdf", sample = sample),
	    "rseqcOutput/GSE110042.geneBodyCoverage.r"

################################################################
## rRNA count
################################################################
rule rRNA_count:
	input:
		"starOutput/{sample}_Aligned.sortedByCoord.out.bam"
	output:
		"rseqcOutput/{sample}.rRNA.txt"
	params:
		split_bam = "/data/home/lingzili/genomeTools/RSeQC-3.0.0/scripts/split_bam.py",
		rRNA_bed = "/data/home/lingzili/genomeTools/RSeQC-3.0.0/genomeBED/mm10_rRNA.bed",
		outprefix = "{sample}" + ".rRNA"
	priority: 40
	message:
		"Count ribosomal RNA on {wildcards.sample}"
	shell:
		'''
		python {params.split_bam} -i {input} -r {params.rRNA_bed} -o {params.outprefix} 1>{output} &&
		rm *.rRNA.*.bam
        '''
################################################################
## Read distribution
################################################################
rule read_distribution:
    input:
        "starOutput/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        "rseqcOutput/{sample}.readDist.txt"
    params:
	    read_distribution = "/data/home/lingzili/genomeTools/RSeQC-3.0.0/scripts/read_distribution.py",
	    RefBED="/data/home/lingzili/genomeTools/RSeQC-3.0.0/genomeBED/mm10_RefSeq.bed"
    priority: 30
    message: 
	    "Running RseQC read distribution on {wildcards.sample}"
    shell:
        "python {params.read_distribution} -i {input} -r {params.RefBED} 1>{output}"

################################################################
## Insert size
################################################################
rule collect_insert_size:
    input:
        "starOutput/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        "rseqcOutput/{sample}.size.pdf"
    message: "Collecting insert size for {wildcards.sample}"
    priority: 20
    shell:
        "java -jar /data/home/lingzili/genomeTools/picard.jar CollectInsertSizeMetrics \
		I={input} \
		O=rseqcOutput/{wildcards.sample}.size.txt \
		H={output} \
		INCLUDE_DUPLICATES=true \
		HISTOGRAM_WIDTH=800 \
		M=0.1"

################################################################
## Gene body coverage
################################################################
rule geneBody_coverage: 
    input: 
	    "starOutput/"
    output:
        "rseqcOutput/GSE110042.geneBodyCoverage.curves.pdf",
        "rseqcOutput/GSE110042.geneBodyCoverage.r"
    params:
	    geneBody_coverage = "/data/home/lingzili/genomeTools/RSeQC-3.0.0/scripts/geneBody_coverage.py",
	    housekeepBED = "/data/home/lingzili/genomeTools/RSeQC-3.0.0/genomeBED/mm10.HouseKeepingGenes.bed",
	    outprefix = "rseqcOutput/GSE110042"
    message: "Gene body coverage"
	shell:
		'''
		python {params.geneBody_coverage} -r {params.housekeepBED} -i {input} -o {params.outprefix} &&
		mv log.txt rseqcOutput/
        '''