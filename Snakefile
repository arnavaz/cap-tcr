#### This Snakemake file is to generate mixcr clones tracking figures from the fastqs###
mport os.path as path

#### This Snakemake file is to generate mixcr clones tracking figures from the fastqs###
configfile: "config.yaml"
basedir = config['basedir']
inputdir = config['inputdir']
aligndir = config['aligndir']
clonetrack = config['clonetrack']
logdir = config['logdir']
pars_log = config['pars_log']
mixcr = config['mixcr_path']
mixcr_qc = config['mixcr_qc'] 
clonesdir = basedir + "/clones/"
assembledir = basedir + "/assemble"
clonetracker = basedir + "/source/clonetrack.R"
print(clonetracker)
###values for running mixcr### these should be set in the config file###

java_align = "java -jar  " + mixcr + "  align -p rna-seq -s hsa -OallowPartialAlignments=true -OvParameters.geneFeatureToAlign=VGeneWithP -r "
assemble1 = "java -jar   " + mixcr + " assemblePartial "
#assembleEx = "java -jar   " + mixcr + " extendAlignments "
assembleEx = "java -jar   " + mixcr + " mergeAlignments "
assemble_f = "java -jar  " + mixcr + " assemble -r"
export = "java -jar  " + mixcr + " exportClones"

IDS1, = glob_wildcards(path.join(inputdir,"{id}_R1.fastq.gz"))
#IDS2, = glob_wildcards(path.join(inputdir,"{id}_R2.fastq.gz"))

print(IDS1)

###  Rules ####

rule all:
    input:
        expand(path.join(aligndir,"alignments_{id}.vdjca"),id=IDS1),
        expand(path.join(assembledir,"rescue1_{id}.vdjca"),id=IDS1),
        expand(path.join(assembledir,"rescue2_{id}.vdjca"),id=IDS1),
        expand(path.join(assembledir,"extended_{id}.vdjca"),id=IDS1),
        expand(path.join(assembledir,"{id}.clns"),id=IDS1),
        expand(path.join(logdir,"log_assemble_{id}.txt"),id=IDS1),
        expand(path.join(clonesdir,"CLONES_TRA_{id}.txt"),id=IDS1),
        expand(path.join(clonesdir,"CLONES_TRB_{id}.txt"),id=IDS1),
        expand(path.join(logdir,"assemble_stats.csv"),id=IDS1),
        #expand(path.join(logdir,"align_stats.csv"),id=IDS1),
        #expand(path.join(basedir,"project_QC.pdf"),id=IDS1),
        expand(path.join(clonesdir,"clonetrack_{param}_{chain}.pdf"),param=config['pattern'],chain=config['chain'])

#rule biopython:
#    input:
#        'conda_env/biopython_installed'
#    output:
#        'conda_env/dependencies'
#
#    conda: 'conda_env/biopython.yml'
#
#    shell:
#        'touch {output}'



rule extract:
     input:
        R1 = path.join(inputdir,"{id}_R1.fastq.gz"),
        R2 = path.join(inputdir,"{id}_R2.fastq.gz")
     output:
        R1 = path.join(basedir,"{id}_barcode_R1.fastq"),
        R2 = path.join(basedir,"{id}_barcode_R2.fastq")
     params:
        outprefix = lambda wildcards, output: output.R1.split('_barcode_')[0],
        #outprefix = lambda wildcards, {input.R1}.split('_barcode_')[0],
        barcodes = config['barcodes'],
        extract_barcode = config['extract_barcodes_path']
     shell:
        "python {params.extract_barcode} --read1 {input.R1} --read2 {input.R2} --outfile {params.outprefix} --blist {params.barcodes}"

rule align:
     input:
        R1 = path.join(basedir,"{id}_barcode_R1.fastq"),
        R2 = path.join(basedir,"{id}_barcode_R2.fastq"),
     output:
        vdjca = path.join(aligndir,"alignments_{id}.vdjca")
     shell:
        "{java_align} {input.R1} {input.R2} {output.vdjca}"

rule assemble1:
     input:
        vdjca = path.join(aligndir,"alignments_{id}.vdjca")
     output:
        vdjca_res1 = path.join(assembledir,"rescue1_{id}.vdjca")
     shell:
        "{assemble1} {input.vdjca} {output.vdjca_res1}"

rule assemble2:
    input:
        vdjca_res1 = path.join(assembledir,"rescue1_{id}.vdjca")
    output:
        vdjca_res2 = path.join(assembledir,"rescue2_{id}.vdjca")
    shell:
        "{assemble1} {input.vdjca_res1} {output.vdjca_res2}"

rule assembleEx:
    input:
        vdjca_res2 = path.join(assembledir,"rescue2_{id}.vdjca")
    output:
        vdjca_ext = path.join(assembledir,"extended_{id}.vdjca")
    shell:
        "{assembleEx} {input.vdjca_res2} {output.vdjca_ext}"

rule assemble:
    input:
        vdjca_ext = path.join(assembledir,"extended_{id}.vdjca")
    output:
        clns = path.join(assembledir,"{id}.clns"),
        log = path.join(logdir,"log_assemble_{id}.txt")
    shell:
        "{assemble_f} {output.log} {input.vdjca_ext} {output.clns}"

rule export_A:
    input:
        cln = path.join(assembledir,"{id}.clns")
    output:
        clones_A = path.join(clonesdir,"CLONES_TRA_{id}.txt")
    shell:
        "{export} --chains TRA {input.cln}  {output.clones_A}"

rule export_B:
    input:
        clns = path.join(assembledir,"{id}.clns")
    output:
        clones_B = path.join(clonesdir,"CLONES_TRB_{id}.txt")
    shell:
        "{export} --chains TRB {input.clns} {output.clones_B}"
rule log_pars:
    output:
        a1 = path.join(logdir,"assemble_stats.csv"),
        #a2 = path.join(logdir,"align_stats.csv")
    shell:
        "python {pars_log} -i  {logdir}  -o {logdir} -pm log_assemble_"

rule run_QC_plotter:
    output:
        qc = path.join(basedir,"project_QC.pdf")
    shell:
        "Rscript {mixcr_qc} {logdir}"

rule clone_tracker:
    input:
        clones = clonesdir
    output:
        #clonetrack = path.join(clonesdir,"clonetrack_.{config['pattern']}.{config['chain']}.pdf")
        clonetrack =  expand(path.join(clonesdir,"clonetrack_{param}_{chain}.pdf"),param=config['pattern'],chain=config['chain'])
    params:
        param1 = config['pattern'],
        type = config['chain']
    shell:
        "Rscript {clonetracker} {input.clones} {params.param1} {params.type}"
