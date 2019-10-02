# Pipeline for trimming, STAR alignment, and finally fastqc for fastq and BAM files
# Command for dry run: snakemake --snakefile star.snakefile -n -p -j 20
# Command to execute: snakemake --snakefile star.snakefile -j 20
# Remember -j will replace the specific threads defined in the rules!

os.system("mkdir -p trimOutput")
os.system("mkdir -p starOutput")
os.system("mkdir -p fastqcOutput")

sample = ["SRR6664591", "SRR6664592", "SRR6664615", "SRR6664616"]

rule all:
    input:
	    "fastqcOutput/multiqc.html",
	    expand("starOutput/{sample}_Aligned.sortedByCoord.out.bam.bai", sample = sample)
	    
################################################################
## Read trimming
################################################################
rule trim:
	input:
		"fastqFile/{sample}_1.fastq.gz", 
		"fastqFile/{sample}_2.fastq.gz"
	output:
		"trimOutput/{sample}_1_val_1.fq.gz",
		"trimOutput/{sample}_2_val_2.fq.gz",
		"trimOutput/{sample}_1.fastq.gz_trimming_report.txt",
		"trimOutput/{sample}_2.fastq.gz_trimming_report.txt"
	message:"Trimming reads on sample {sample}"
	threads: 10
	shell:
		"trim_galore --cores {threads} --paired -o trimOutput {input}"

################################################################
## STAR alignment
################################################################
rule align_STAR:
	input:
		read1 = "trimOutput/{sample}_1_val_1.fq.gz",
		read2 = "trimOutput/{sample}_2_val_2.fq.gz",
		STARindex = "/data/Genomes/mouse/mm10/star/index101bp/"
	output:
		"starOutput/{sample}_Aligned.sortedByCoord.out.bam",
		"starOutput/{sample}_Log.final.out"
	params:
		outprefix = "{sample}" + "_"
	message: "Running STAR Alignment on {sample}"
	threads: 10
	shell:
		"/data/home/lingzili/.conda/envs/lingzili/bin/STAR --runThreadN {threads} \
		--genomeDir {input.STARindex} \
		--readFilesIn {input.read1} {input.read2} \
		--readFilesCommand zcat \
		--outFileNamePrefix ./starOutput/{params.outprefix} \
		--outSAMtype BAM SortedByCoordinate"

rule index_bam:
	input:
		"starOutput/{sample}_Aligned.sortedByCoord.out.bam"
	output:
		"starOutput/{sample}_Aligned.sortedByCoord.out.bam.bai"
	message: "Indexing BAM files of {sample}"
	threads: 10
	shell:
		"samtools index -@ {threads} {input} {output}"

################################################################
## FastQC and MultiQC
################################################################
rule fastqc: 
    input: 
        "trimOutput/{sample}_1_val_1.fq.gz", 
        "trimOutput/{sample}_2_val_2.fq.gz",
        "starOutput/{sample}_Aligned.sortedByCoord.out.bam"
    output: 
        "fastqcOutput/{sample}_1_val_1_fastqc.zip",
        "fastqcOutput/{sample}_2_val_2_fastqc.zip",
        "fastqcOutput/{sample}_Aligned.sortedByCoord.out_fastqc.zip"
    message: "FastQC"
    threads: 10
    shell:
        "fastqc -t {threads} -o fastqcOutput {input}"

rule multiqc: 
    input: 
	    expand(["trimOutput/{sample}_1.fastq.gz_trimming_report.txt", "trimOutput/{sample}_2.fastq.gz_trimming_report.txt",
	    	"fastqcOutput/{sample}_1_val_1_fastqc.zip", "fastqcOutput/{sample}_2_val_2_fastqc.zip",
	    	"fastqcOutput/{sample}_Aligned.sortedByCoord.out_fastqc.zip", "starOutput/{sample}_Log.final.out"], sample = sample)
    output:
        "fastqcOutput/multiqc.html"
    message: "MultiQC"
	shell:
		"/usr/local/bin/multiqc -o fastqcOutput {input}"
