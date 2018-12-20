shell.prefix("source activate mageck_vispr")
DIR_RAW = config["raw_data_dir"]
SAMPLE = glob_wildcards(DIR_RAW + "/{sample}_1.fastq").sample

DIR_FQ = DIR_RAW
TRIM_GALORE = "s1_trim_galore/"
TRIM_GALORE_FQC = "s1_trim_galore_fastqc/"
LINKER_RM_DIR="s2_linker_rm/"
FASTP_SPLIT="s2_fastp_split/"
BOWTIE2_ALN="s3_bowtie2_aln/"
BOWTIE2_ALN_MER="s3_bowtie2_aln_mer/"
HICBuildMatrix = "s4_hicBuildMatrix/"

FILE_NUM = 10
if "fastp_split" in config and "file_num" in  config["fastp_split"]:
    FILE_NUM = config["fastp_split"]["file_num"]

PREFIX = []
for i in range(1, FILE_NUM + 1):
    PREFIX.append('{num:04d}'.format(num=i))

print(PREFIX)
#https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#local-rules
localrules: all, trim_galore, fastqc_trim_galore, linker_rm_ng, bowtie2_aln_merge

rule all:
    input:
        #expand(DIR_FQ + "{sample}_{group}.fastq", sample=SAMPLE, group=[1,2]),
        expand(TRIM_GALORE_FQC + "{sample}_1_trimmed_fastqc.zip", sample=SAMPLE),
        expand(LINKER_RM_DIR + "{sample}_{group}.fq", sample=SAMPLE, group=["left","right"]),
        expand(FASTP_SPLIT + "fastp_split.{sample}.done", sample=SAMPLE),
        expand(BOWTIE2_ALN + "{prefix}.{sample}.{lr}.bam", prefix =PREFIX, sample=SAMPLE, lr=["right", "left"]),
        expand(BOWTIE2_ALN_MER + "{sample}.{lr}.bam", sample=SAMPLE, lr=["right", "left"]),
        expand(HICBuildMatrix + "{sample}.hicBuildMatrix.bam", sample=SAMPLE),
        expand(HICBuildMatrix + "{sample}.plotMatrix.png", sample=SAMPLE),

rule trim_galore:
    input:
        read1=DIR_FQ + "{sample}_1.fastq",
        #read2=DIR_FQ + "{sample}_2.fastq",
    #group: "read_processing"
    params: TRIM_GALORE
    output:
        out_fq1=TRIM_GALORE + "{sample}_1_trimmed.fq",
        #out_fq2=TRIM_GALORE + "{sample}_2_trimmed.fq"
    shell:
        '''
trim_galore --length 75 {input.read1} --output_dir {params}
'''

rule fastqc_trim_galore:
    input:
        read1=TRIM_GALORE + "{sample}_1_trimmed.fq",
        #read2=TRIM_GALORE + "{sample}_2_trimmed.fq",
    #group: "read_processing"
    params: TRIM_GALORE_FQC
    output:
        out_fq1=TRIM_GALORE_FQC + "{sample}_1_trimmed_fastqc.zip",
        #out_fq2=TRIM_GALORE_FQC + "{sample}_2_trimmed.fq_fastqc.zip"
    shell:
        '''
fastqc {input.read1} -o {params}
'''

rule linker_rm_ng:
    input:
        read1=TRIM_GALORE + "{sample}_1_trimmed.fq",
        #read2=TRIM_GALORE + "{sample}_2_trimmed.fq",
    #group: "read_processing"
    params:
        script=config["scripts"]["linker_rm"],
        adapter=config["adapter"],
        prefix=LINKER_RM_DIR + "{sample}"
    output:
        out1_fq1=LINKER_RM_DIR+"{sample}_left.fq",
        out1_fq2=LINKER_RM_DIR+"{sample}_right.fq",
        #out2_fq1=LINKER_RM_DIR+"{sample}_read2_R1.fq",
        #out2_fq2=LINKER_RM_DIR+"{sample}_read2_R2.fq",
        #read_left=LINKER_RM_DIR+"{sample}_left.fq",
        #read_right=LINKER_RM_DIR+"{sample}_right.fq",
    shell:
        '''
python {params.script} --read {input.read1} -l {params.adapter} -c 70 -p {params.prefix}
'''

localrules: fastp_split
rule fastp_split:
    input:
        read_left=LINKER_RM_DIR+"{sample}_left.fq",
        read_right=LINKER_RM_DIR+"{sample}_right.fq",
    #group: "read_processing"
    params: 
        file_num = FILE_NUM,
        suffix_left="{sample}.left.fq",
        suffix_right="{sample}.right.fq",
        out_dir=FASTP_SPLIT
    output:
        touch_done=FASTP_SPLIT + "fastp_split.{sample}.done"
        #fq=expand(FASTP_SPLIT + "{prefix}.{sample1}_{lr}.fq", prefix =PREFIX, sample1=SAMPLE, lr=["right", "left"])
    shell:
        '''
fastp -i {input.read_left} -o {params.suffix_left} --split {params.file_num}
mv *{params.suffix_left} {params.out_dir}
fastp -i {input.read_right} -o {params.suffix_right} --split {params.file_num}
mv *{params.suffix_right} {params.out_dir}
touch {output.touch_done}
'''

CPU_NUM = 1
if "bowtie2_aln" in config and "cpu_num" in config["bowtie2_aln"]:
    #print(config["bowtie2_aln"]["cpu_num"])
    CPU_NUM = config["bowtie2_aln"]["cpu_num"]
print("CPU_NUM:" +str(CPU_NUM))

rule bowtie2_aln_split:
    input:
        touch_done=FASTP_SPLIT + "fastp_split.{sample}.done" 
        #read_left=FASTP_SPLIT+"{prefix}.{sample}.left.fq",
        #read_right=FASTP_SPLIT+"{prefix}.{sample}.right.fq"
    #group: "read_processing"
    log:
        left=BOWTIE2_ALN + "{prefix}.{sample}.left.log",
        right=BOWTIE2_ALN + "{prefix}.{sample}.right.log", 
    params:
        read_left=FASTP_SPLIT+"{prefix}.{sample}.left.fq",
        read_right=FASTP_SPLIT+"{prefix}.{sample}.right.fq", 
        cpu_num=CPU_NUM,
        reference=config["bowtie2_aln"]["reference"],
    output:
        bam_left=BOWTIE2_ALN + "{prefix}.{sample}.left.bam",
        bam_right=BOWTIE2_ALN + "{prefix}.{sample}.right.bam"
    shell:
        '''
bowtie2 -p {params.cpu_num} --reorder -D 15 -R 2 -N 0 -L 8 -i S,1,0.75 -x {params.reference} -U {params.read_left} 2> {log.left} | samtools view -Shb - -o {output.bam_left}
bowtie2 -p {params.cpu_num} --reorder -D 15 -R 2 -N 0 -L 8 -i S,1,0.75 -x {params.reference} -U {params.read_right} 2> {log.right} | samtools view -Shb - -o {output.bam_right}
'''
rule bowtie2_aln_merge:
    input:
        bam_left=expand(BOWTIE2_ALN + "{prefix}.{sample}.left.bam", prefix=PREFIX, sample=SAMPLE), 
        bam_right=expand(BOWTIE2_ALN + "{prefix}.{sample}.left.bam", prefix=PREFIX, sample=SAMPLE),
    #group: "samtools_merge"
    params: 
        cpu_num=CPU_NUM
    output:
        bam_mer_left=BOWTIE2_ALN_MER + "{sample}.left.bam",
        bam_mer_right=BOWTIE2_ALN_MER + "{sample}.right.bam"
    shell:
        '''
samtools merge -@ {params.cpu_num} {output.bam_mer_left} {input.bam_left} 
samtools merge -@ {params.cpu_num} {output.bam_mer_right} {input.bam_right}
'''


localrules: hicBuildMatrix

rule hicBuildMatrix:
    input:
        bam_left=BOWTIE2_ALN_MER + "{sample}.left.bam",
        bam_right=BOWTIE2_ALN_MER + "{sample}.right.bam",
    #group: "hicexplorer"
    params: 
        QCfolder=HICBuildMatrix + "QCfolder_{sample}",
        binSize = config["hicBuildMatrix"]["binSize"],
        #bufferSize=config["hicBuildMatrix"]["bufferSize"]
    log: 
        buildMatrix=HICBuildMatrix + "{sample}_hicBuildMatrix.log",
        plotMatrix=HICBuildMatrix + "{sample}_res_hicPlotMatrix.23chr.bin.log"
    threads: 14
    output:
        bam=HICBuildMatrix + "{sample}.hicBuildMatrix.bam",
        h5=HICBuildMatrix + "{sample}.hicBuildMatrix.h5",
        h5_mer=HICBuildMatrix + "{sample}.hicBuildMatrix.merge.h5",
        plotMatrix=HICBuildMatrix + "{sample}.plotMatrix.png",
    shell:
        '''
time hicBuildMatrix --samFiles {input.bam_left} {input.bam_right} --binSize {params.binSize} --restrictionSequence AATT --threads {threads}  --outBam {output.bam} -o {output.h5} --QCfolder {params.QCfolder} 2>{log.buildMatrix}
time hicMergeMatrixBins -m {output.h5} -o {output.h5_mer} -nb 100
time hicPlotMatrix --matrix {output.h5_mer} --outFileName {output.plotMatrix} --log1p --dpi 900 --colorMap Reds --chromosomeOrder chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY --title chr23 2> {log.plotMatrix}
'''
