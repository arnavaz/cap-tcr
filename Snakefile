#### This Snakemake file is to generate mixcr clones tracking figures from the fastqs###
configfile: "config.yaml"
basedir = config['basedir']
inputdir = config['inputdir']
aligndir=config['aligndir']
clonetrack=config['clonetrack']
logdir=config['logdir']
pars_log=config['pars_log_path']
mixcr=config['mixcr_path']
clonesdir=basedir + "/clones"
assembledir= basedir + "/assemble"
###values for running mixcr### these should be set in the config file###

java_align = "java -jar  " +mixcr+ "  align -p rna-seq -s hsa -OallowPartialAlignments=true -OvParameters.geneFeatureToAlign=VGeneWithP -r "
assemble1 = "java -jar   "+ mixcr+ " assemblePartial " 
assembleEx="java -jar   "+ mixcr + " extendAlignments "
assemble_f="java -jar  "+mixcr + " assemble -r"
export="java -jar  "+ mixcr + " exportClones"

IDS1, = glob_wildcards(path.join(inputdir,"{id}_R1.fastq.gz"))
IDS2, = glob_wildcards(path.join(inputdir,"{id}_R2.fastq.gz"))
#print("test")

###  Rules ####

rule all:
            input: 
                   #expand("clonesdir/clonetrack_{id}.pdf",id=IDS1)
		   #expand('clonesdir/clonetrack_{id}.{param}.pdf'.format(param=config['pattern']),id=IDS1)
                   expand(path.join(clonesdir,"clonetrack_{id}"), param=config['pattern'], id=IDS1)
		
#rule biopython:
#        output:
#               'conda_env/biopython_installed'
#        conda: 'conda_env/biopython.yml'
#        shell:
#               'touch {output}'

rule extract:
	    input: 
                   R1 = exapnd(path.join(inputdir,"{id}_R1.fastq.gz"),id=IDS1), 
	           R2 = expand(path.join(inputdir,"{id}_R2.fastq.gz"),id=IDS2)
            output:
		    path.join(basedir,"{id}_barcode_R1.fastq"),
	            path.join(basedir,"{id}_barcode_R2.fastq")

            params:
                   #outprefix = lambda wildcards,  {input.R1}.split('_barcode_')[0],
                   outprefix="JAX",
	           barcodes = config['barcodes'],
	           extract_barcode = config['extract_barcodes_path']		
            shell:
                  "python {params.extract_barcodes} --read1 {input.R1} --read2 {input.R2} --outfile {params.outprefix} {params.barcodes}"
 

rule align:

	    input:
	      R1 = expand("basedir/{id}_barcode_R1.fastq", id=IDS1),
	      R2 = expand("basedir/{id}_barcode_R1.fastq", id=IDS2)
	
	    output: 
		   vdjca = "aligndir/alignments_{id}.vdjca"
	    shell:
		   "{java_align} logdir/{id}.txt {input.R1} {input.R2} {output.vdjca}"
		   


rule assemble1:

	    input:
	          vdjca = "aligndir/alignments_{id}.vdjca"

	    output:
	          vdjca_res1= "assembledir/rescue1_{id}.vdjca"
		
	    shell:
                   "{assembl1} {input.vdjca} {output.vdjca_res1}"



rule assemble2:

	    input:
	          vdjca_res1 = "rescue1_{id}.vdjca"

	    output:
	          vdjca_res2="assembledir/rescue2_{id}.vdjca"

	    shell:
                  "{assemble} {input.vdjca_res1} {output.vdjca_res2}"

rule assembleEx:

	    input:
		vdjca_res2="assembledir/rescue2_{id}.vdjca"

	    output:
               vdjca_ext="assembledir/extended_{id}.vdjca"
	    shell:
	          "{assembleEx} {input.vdjca_res2} {output.vdjca_ext}"


rule assemble:

	    input:
		vdjca_ext="assembledir/extended_{id}.vdjca"
	    output:
	        clns="assembledir/{id}.clns",
		log="logdir/log_assemble_id}.txt"			    

	    shell:
	          "{assemble_f} -r logdir/log_assemble_{id}.txt assembledir/extended_{id}.vdjca {output.clns}"
	          


rule export_A:

	    input:
		   cln="assembledir/{id}.clns"
		
	    output:
          	   clones_A="clonesdir/CLONES_TRA_{id}.txt"	
		
	
	    shell:
		 "export exportClones --chains TRA {input.cln}  {output.clones_A}"
rule export_B:
	
	    input:
		   clns="assembledir/{id}.clns"
	    output:
          	   clones_B="clonesdir/CLONES_TRB_{id}.txt"		
	
            shell:
		 "export  exportClones --chains TRB assembledir/{input.clns} {output.clones_B}"

rule log_pars:
	
	    output: 
		a1="logdir/assemble_stats.csv",
		a2="logdir/align_stats.csv"
	    shell:
		 "python {parse_log_path} -i  logdir  -o logdir"
		
rule run_QC_plotter:
	
	    input:
		a1="logdir/assemble_stats.csv",
		a2="logdir/align_stats.csv"
		
	    output:
		qc="project_QC.pdf"
		
	    shell:
		  "Rscript {mixcr_qc}  {input.a1} {input.a2}"
rule  clone_tracker:
	
	    input:
		clones='clonesdir/CLONES_{chain}_{id}.txt'.format(chain=config['chain'])
		
	    output:
		clonetrack='clonesdir/clonetrack_{id}.{param}.pdf'.format(param=config['pattern'])
            params:
                   param1=config['pattern']
		   type=config['chain']		
	    shell:
	          "Rscript {clonetracker} clonesdir {params.param1}  {params.type}"
		   	
