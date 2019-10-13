# Command for dry run: snakemake --snakefile mm10.ATAC.snakefile -n -p
# Command for dag: snakemake --snakefile mm10.ATAC.snakefile --dag | dot -Tpdf > dag.pdf
# Command for execute: snakemake --snakefile mm10.ATAC.snakefile -j 20 -w 30

# Specify fastq file location
vds2Dir = "/data/home/lingzili/ATAC_analysis/19032019/fastqFile"

# Extract sample ID from pair-end fastq files
import os

extractList = []

for fileName in os.listdir(vds2Dir):
	extractList.append(fileName.split('_')[0])

# Get unique sample ID
ID = list(set(extractList))

# Create directories
os.system("mkdir -p trimOutput")
os.system("mkdir -p Bowtie2Output")
os.system("mkdir -p peakOutput")
os.system("mkdir -p fastqcOutput")
os.system("mkdir -p ATACqcOutput")

rule final_results:
    input:
	    expand("peakOutput/{sample}.narrowPeak", sample = ID),
	    expand("peakOutput/{sample}.bedgraph", sample = ID),
	    expand("ATACqcOutput/{sample}.phantom.pdf", sample = ID),
	    expand("ATACqcOutput/{sample}.lib.metrics.txt", sample = ID),
	    "fastqcOutput/multiqc_report.html"

################################################################
## Read trimming
################################################################
rule trim:
	input:
		vds2Dir + "/{sample}_R1.fastq.gz", 
		vds2Dir + "/{sample}_R2.fastq.gz"
	output:
		"trimOutput/{sample}_R1_val_1.fq.gz",
		"trimOutput/{sample}_R2_val_2.fq.gz",
		"trimOutput/{sample}_R1.fastq.gz_trimming_report.txt",
		"trimOutput/{sample}_R2.fastq.gz_trimming_report.txt"
	message:"Trimming reads on sample {wildcards.sample}"
	threads: 6
	shell:
		"trim_galore --cores {threads} --paired -o trimOutput {input}"

################################################################
## Bowtie2 alignment
################################################################
rule align_Bowtie2:
	input:
		read1 = "trimOutput/{sample}_R1_val_1.fq.gz", 
		read2 = "trimOutput/{sample}_R2_val_2.fq.gz"
	output:
		"Bowtie2Output/{sample}.bam"
	params:
	    Bowtie2index = "/data/home/lingzili/mm10_genome/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome"
	log:
	    "Bowtie2Output/{sample}.bam.log"
	message: "Running Bowtie2 alignment on {wildcards.sample} and sort by coordinate"
	threads: 10
	shell:
		"bowtie2 -p {threads} --very-sensitive -X 2000 -x {params.Bowtie2index} -1 {input.read1} -2 {input.read2} 2>{log} | samtools view -bSh -@ {threads} \
		| samtools sort -m 1G -@ {threads} -O bam -o {output} && samtools index -@ {threads} {output}"

rule flagstat_bam:
	input:
		"Bowtie2Output/{sample}.bam"
	output:
		"Bowtie2Output/{sample}.bam.flagstat.txt"
	message: "Flagstat on {wildcards.sample}"
	threads: 10
	shell:
		"samtools flagstat -@ {threads} {input} 1>{output}"

################################################################
## Read distribution
################################################################
rule read_distribution:
    input:
        "Bowtie2Output/{sample}.bam"
    output:
        "ATACqcOutput/{sample}.readDist.txt"
    params:
	    read_distribution = "/data/home/lingzili/genomeTools/RSeQC-3.0.0/scripts/read_distribution.py",
	    RefBED = "/data/home/lingzili/genomeTools/RSeQC-3.0.0/genomeBED/mm10_RefSeq.bed"
    message: 
	    "Running read distribution on {wildcards.sample}"
    shell:
        "python {params.read_distribution} -i {input} -r {params.RefBED} 1>{output}"

################################################################
## Library complexity
################################################################
rule library_complexity:
    input:
        "Bowtie2Output/{sample}.bam"
    output:
        "ATACqcOutput/{sample}.lib.metrics.txt"
    message: "Library complexity for {wildcards.sample}"
    shell:
        "java -jar /data/home/lingzili/genomeTools/picard.jar EstimateLibraryComplexity I={input} O={output}"

################################################################
## Post-alignment filtering
################################################################
# Locate and tag duplicates in BAM file
rule mark_duplicates:
	input:
		"Bowtie2Output/{sample}.bam"
	output:
		bam = "Bowtie2Output/{sample}.markDup.bam",
		dupMetrics = "Bowtie2Output/{sample}.dupMetrics.txt"
	message: 
	    "Mark duplicates on {wildcards.sample} BAM file"
	shell:
		"java -jar /data/home/lingzili/genomeTools/picard.jar MarkDuplicates \
		I={input} O={output.bam} M={output.dupMetrics} ASSUME_SORTED=true REMOVE_DUPLICATES=false QUIET=true VERBOSITY=ERROR"

# -F 1804: Remove read unmapped(0x4), mate unmapped(0x8), not primary alignment(0x100), read fails quality checks(0x200) and PCR duplicate (0x400)
# -f 0x2: Only output alignments with read mapped in proper pair(0x2)
# -q 30: Skip alignments with MAPQ<30
# grep -v chrM: Print reads that do not contain chrM
rule filter_reads:
	input:
		"Bowtie2Output/{sample}.markDup.bam"
	output:
		"Bowtie2Output/{sample}.flt.sortByName.bam"
	message: "Filter reads on {wildcards.sample} and sort by name"
	threads: 10
	shell:
		"samtools view -@ {threads} -F 1804 -f 0x2 -q 30 -h {input} | grep -v chrM | samtools sort -@ {threads} -n -m 1G -O bam -o {output}"

################################################################
## Peak calling by Genrich
################################################################
# ATAC-seq mode (-j)
# Remove chrY (-e chrY)
# Keep unpaired alignments (-y)
# Remove PCR duplicates (-r)
# Print status updates/counts to stderr (-v)
rule peak_calling:
    input:
	    "Bowtie2Output/{sample}.flt.sortByName.bam"
    output:
        np = "peakOutput/{sample}.narrowPeak",
        log = "peakOutput/{sample}.peak.log"
    log:
	    "peakOutput/{sample}.Genrich.log"
    message: "Peak calling on {wildcards.sample}"
    shell:
	    "Genrich -t {input} -j -e chrY -y -r -v -o {output.np} -f {output.log} 2>{log}"

################################################################
## bamCoverage
################################################################
rule bam_coverage:
	input:
		"Bowtie2Output/{sample}.bam"
	output:
		"peakOutput/{sample}.bedgraph"
	message: "Create bedGraph of {wildcards.sample}"
	threads: 10
	shell:
		"bamCoverage -p {threads} -b {input} -of bedgraph -o {output} \
		-bs 10 --extendReads --ignoreDuplicates --normalizeUsing RPKM --samFlagExclude 1804"

################################################################
## Phantom peak
################################################################
rule phantom_peak:
    input:
	    "Bowtie2Output/{sample}.flt.sortByName.bam"
    output:
        phantomTxt = "ATACqcOutput/{sample}.phantom.txt",
        plot = "ATACqcOutput/{sample}.phantom.pdf"
    log:
		"ATACqcOutput/{sample}.phantom.log"
    params:
	    run_spp = "/data/home/lingzili/genomeTools/phantompeakqualtools/run_spp.R"
    threads: 10
    message: "Phantom peaks for {wildcards.sample}"
    shell: """ 
        /data/home/lingzili/.conda/envs/lingzili/bin/Rscript {params.run_spp} -p={threads} -c={input} -odir=ATACqcOutput -savp={output.plot} -out={output.phantomTxt} 2>{log}
		"""

################################################################
## FastQC and MultiQC
################################################################
rule fastqc: 
    input: 
        "trimOutput/{sample}_R1_val_1.fq.gz", 
        "trimOutput/{sample}_R2_val_2.fq.gz"
    output: 
        "fastqcOutput/{sample}_R1_val_1_fastqc.zip",
        "fastqcOutput/{sample}_R2_val_2_fastqc.zip"
    message: "FastQC"
    threads: 10
    shell:
        "fastqc -t {threads} --quiet -o fastqcOutput {input}"

rule multiqc: 
    input: 
	    expand(["trimOutput/{sample}_R1.fastq.gz_trimming_report.txt", "trimOutput/{sample}_R2.fastq.gz_trimming_report.txt", 
	    	"fastqcOutput/{sample}_R1_val_1_fastqc.zip", "fastqcOutput/{sample}_R2_val_2_fastqc.zip", 
	    	"Bowtie2Output/{sample}.dupMetrics.txt", "Bowtie2Output/{sample}.bam.log", "Bowtie2Output/{sample}.bam.flagstat.txt", 
	    	"ATACqcOutput/{sample}.readDist.txt"], sample = ID)
    output:
        "fastqcOutput/multiqc_report.html"
    message: "MultiQC"
	shell:
		"/usr/local/bin/multiqc -o fastqcOutput {input}"
