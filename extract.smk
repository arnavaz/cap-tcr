
import os
import sys
configfile: "config.yaml"
basedir = config['basedir']
inputdir = config['inputdir']
IDS1, = glob_wildcards("inputdir/{id}_R1.fastq.gz")
IDS2, = glob_wildcards("inputdir/{id}_R2.fastq.gz")


rule all:
    input:  expand("basedir/{id}_barcode_R1.fastq", id=IDS1),
	    expand("basedir/{id}_barcode_R2.fastq", id=IDS2)


rule extract:

        conda: 'conda_env/biopython.yml'
	input: 
              R1=expand("inputdir/{id}_R1.fastq.gz", id=IDS1), 
	      R2=expand("inputdir/{id}_R2.fastq.gz",id=IDS2)
        output:
	     R1 = expand("basedir/{id}_barcode_R1.fastq", id=IDS1),
	     R2 = expand("basedir/{id}_barcode_R1.fastq", id=IDS2)

       params:
            outprefix = lambda wildcards, output: output.R1.split('_barcode_')[0]
            barcodes = config['barcodes']
       shell:
        '''
        python {extract_barcodes}
            --read1 {{input.R1}}
            --read2 {{input.R2}}
            --outfile {{params.outprefix}}
            {{params.barcodes}}
        '''.format(extract_barcodes = config['extract_barcodes_path'])

