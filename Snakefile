"""
mTORC3 RNA-Seq Snakemake pipeline
V. Keith Hughitt
Oct 27, 2020

Basic pipeline to perform clean and analyze human RNA-Seq samples relating to Joe's
mTORC3 project.

Since this pipeline is being used as a teaching opportunity, comments will be included
throughout the file in order to improve clarity.

TODO AT NEXT MEETING:

1. Conda env review (activate)
2. Walk through snakefile changes
3. Snakemake workflow dev process (start simple)
    - Best thing is to first run each step manually outside of Snakemake to make sure
      it works as expected.

For more information, see:

- https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html

"""
import os

#
# configuration
#
data_dir = "/data/raw/mtorc3/fastq"
out_dir = '/data/proj/mtorc3'

#
# The first rule in a pipeline is often one which sites at the top of the workflow graph
# and thus depends on many other rules / comes last.
# It is used to help define what the final expected outputs are, which in turn
# determines the expected inputs/outputs for other dependent rules in the pipeline.
#

# 
# There are different ways one can specific the expected inputs, but one relatively
# simple approach I like to use is to use python to get a list of input files, parse out
# the sample ids (everything but the file extension in the filenames), and then feed
# this into a snakemake "expand()" call.
# See the Snakemake docs on expand() for more information.
#

# Sample ids - note that normally we might parse these out from the file names using
# python. for simplicity here, however, they are simply hard-coded below.
samples = ['Rap_d1_starved_2', 'Rap_d1_starved_4', 'Rap_d1_nutri_2', 'Rap_d1_starved_1',
           'Rap_d1_starved_3', 'Rap_d1_nutri_3', 'Ric_d5_nutri_2', 'Rap_d1_nutri_4',
           'Ric_d5_nutri_1', 'Rap_d1_nutri_1']

# wildcards
cell_lines = ['v_d4', 'm7_d5', 'Rap_d1', 'Ric_d5']
nutrients = ['nutri', 'starved']
replicates = ['1', '2', '3', '4']

#
# Each paired-end RNA-Seq sample is associated with two fastq files: "R1" and "R2".
#
rule all:
    input:
        expand(os.path.join(out_dir, 'hisat2', '{cell}_{nutrient}_{replicate}_all_trimmed.sam'),
                  cell=cell_lines, nutrient=nutrients, replicate=replicates)

rule samtools_bam:
    input:
        os.path.join(out_dir, 'hisat2', "{cell}_{nutrient}_{replicate}_all_trimmed.sam")
    output:
        os.path.join(out_dir, 'hisat2', "{cell}_{nutrient}_{replicate}_all_trimmed.bam")
    shell:
        """
        samtools view -S -b {input} > {output}
        """

#
# HISAT2 is a modern RNA-Seq read mapper that is much faster and requires much less
# memory than earlier tools such as tophat and bwa.
#
# Like many other read alignment tools, it requires a reference genome "index", which
# provides useful information about the human reference genome that we will use that is
# independent of the reads we want to map, and thus, can be pre-computed once and
# re-used each time we wish to use HISAT2 for a difference project.
#
# In order to save users to time and effort of computing this index themselves, HISAT2
# provides a link you can use to download indices for the GRCh38 referance genome.
#
# In the rule below, these are assumed to have been already downloaded and places in a
# folder "/data/ref/hisat2/grch38/".
#
# For more information on HISAT2 and what the different parameters below mean, check out
# the manual at: https://daehwankimlab.github.io/hisat2/manual/
#
rule hisat2:
    input:
        r1=os.path.join(out_dir, "cutadapt", "{cell}_{nutrient}_{replicate}_R1_all_trimmed.fastq.gz"),
        r2=os.path.join(out_dir, "cutadapt", "{cell}_{nutrient}_{replicate}_R2_all_trimmed.fastq.gz")
    output:
        os.path.join(out_dir, 'hisat2', "{cell}_{nutrient}_{replicate}_all_trimmed.sam")
    shell:
        """
        hisat2 \
           -x /data/ref/hisat2/grch38/genome \
           -1 {input.r1} \
           -2 {input.r2} \
           --no-unal \
           -S {output}
        """

#
# After we have trimmed our reads, it's always a good idea to run FastQC again on those
# trimmed reads to assess the impact of trimming and to make sure that the quality is
# sufficient before moving onto mapping.
#
rule fastqc_trimmed:
    input:
        r1=os.path.join(out_dir, "cutadapt", "{cell}_{nutrient}_{replicate}_R1_all_trimmed.fastq.gz"),
        r2=os.path.join(out_dir, "cutadapt", "{cell}_{nutrient}_{replicate}_R2_all_trimmed.fastq.gz")
    output:
        os.path.join(out_dir, 'fastqc', 'trimmed', '{cell}_{nutrient}_{replicate}_R1_all_trimmed_fastqc.html'),
        os.path.join(out_dir, 'fastqc', 'trimmed', '{cell}_{nutrient}_{replicate}_R1_all_trimmed_fastqc.zip'),
        os.path.join(out_dir, 'fastqc', 'trimmed', '{cell}_{nutrient}_{replicate}_R2_all_trimmed_fastqc.html'),
        os.path.join(out_dir, 'fastqc', 'trimmed', '{cell}_{nutrient}_{replicate}_R2_all_trimmed_fastqc.zip')
    shell:
        """
        outdir=$(dirname {output[0]})

        fastqc \
            --noextract \
            --quiet \
            --outdir $outdir \
            {input.r1} {input.r2}
        """

#
# After performing an initial check with FastQC, we are ready to perform some basic
# trimming and adaptor removal. It is common for RNA-Seq reads to include some portions
# of the adaptors used in the sample generation. Removing these sequences, which are
# unrelated to the organism's actual genome sequence, will result in a higher percentage
# of reads mapping to the genome, and thus, more information for us to use in our
# downstream analyses. Similarly, the quality of basepair predictions varies across each
# RNA-Seq read, and there is a well-known bias where basepairs at the ends of the reads
# tend be of lower quality (i.e. the fastq file is more likely to report the wrong
# basepair at those positions). This can also lead to a failure to mapping, and as such,
# it is common practice to "trim" low-quality basepairs from the ends of our raw reads.
# cutadapt, used below, helps us with both of these steps, taking our raw reads as input
# and returning a set of "trimmed" reads as output.
#
# For more information, see: https://cutadapt.readthedocs.io/en/stable/
#
rule cutadapt:
    input:
        r1=os.path.join(data_dir, "{cell}_{nutrient}_{replicate}_R1_all.fastq.gz"),
        r2=os.path.join(data_dir, "{cell}_{nutrient}_{replicate}_R2_all.fastq.gz")
    output:
        r1=os.path.join(out_dir, "cutadapt", "{cell}_{nutrient}_{replicate}_R1_all_trimmed.fastq.gz"),
        r2=os.path.join(out_dir, "cutadapt", "{cell}_{nutrient}_{replicate}_R2_all_trimmed.fastq.gz")
    shell:
        """
        cutadapt \
            -o {output.r1} \
            -p {output.r2} \
            -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            -q 20 \
            --minimum-length 25 \
            {input.r1} {input.r2}
        """

#
# The first first step in our pipeline will be to run FastQC on our _raw_ RNA-Seq reads.
# This is useful to get a sense for the quality of the samples, and to check for
# potential problems that we may want to deal with before proceeding to later steps such
# as read mapping and quantification.
#
# For more information on how to run FastQC, and how to interpret its output, see:
# https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
#
rule fastqc_raw:
    input:
        r1=os.path.join(data_dir, "{cell}_{nutrient}_{replicate}_R1_all.fastq.gz"),
        r2=os.path.join(data_dir, "{cell}_{nutrient}_{replicate}_R2_all.fastq.gz")
    output:
        os.path.join(out_dir, 'fastqc', 'untrimmed', '{cell}_{nutrient}_{replicate}_R1_all_fastqc.html'),
        os.path.join(out_dir, 'fastqc', 'untrimmed', '{cell}_{nutrient}_{replicate}_R1_all_fastqc.zip'),
        os.path.join(out_dir, 'fastqc', 'untrimmed', '{cell}_{nutrient}_{replicate}_R2_all_fastqc.html'),
        os.path.join(out_dir, 'fastqc', 'untrimmed', '{cell}_{nutrient}_{replicate}_R2_all_fastqc.zip')
    shell:
        """
        outdir=$(dirname {output[0]})

        fastqc \
            --noextract \
            --quiet \
            --outdir $outdir \
            {input.r1} {input.r2}
        """
