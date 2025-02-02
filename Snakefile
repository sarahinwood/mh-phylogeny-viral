#!/usr/bin/env python3

blast_container = 'docker://ncbi/blast:2.12.0'
raxmlng_container = 'docker://evolbioinfo/raxml-ng:v1.0.2'

#########
# RULES #
#########

rule target:
    input:
        'output/dafv_blastp/dafv_blastp.outfmt6',
        'output/raxml/raxml_all/Blosum62.raxml.bestTree'

######################
## concat phylogeny ##
######################

##draw tree using iTOL and edit in inkscape

##grep "Final LogLikelihood:" *.raxml.log - which model has highest likelihood?
##Isn't a metric to compare between different trees, but between same tree with different settings
rule raxml_all:
    input:
        msa = 'output/raxml/check/check.raxml.reduced.phy'
    output:
        'output/raxml/raxml_all/Blosum62.raxml.bestTree'
    params:
        wd = 'output/raxml/raxml_all',
        prefix = 'output/raxml/raxml_all/Blosum62'
    singularity:
        raxmlng_container
    log:
        'output/logs/raxml_all.log'
    threads:
        20
    shell:
        'raxml-ng '
        '--all '
        '--msa {input.msa} '
        '--model Blosum62 ' ##blosum62 gives highest likelihood, after testing each model
        '--bs-trees 1000 ' #1000 bootstrap replicates
        '--prefix {params.prefix} '
        '--threads {threads} '
        '--seed 11280 '
        '2> {log}'

rule raxml_check:
    input:
        'output/trimal/trimal.fas'
    output:
        'output/raxml/check'
    params:
        wd = 'output/raxml/check',
        prefix = 'output/raxml/check/check'
    singularity:
        raxmlng_container
    threads:
        20
    shell:
        'mkdir {params.wd} || exit 1 ;'
        'raxml-ng '
        '--check '
        '--msa {input} '
        '--model LG ' ##LbFV used this model - need to check though
        '--prefix {params.prefix}'

##trim alignments
rule trimal_concat:
    input:
        'output/muscle/FcC_supermatrix.fas'
    output:
        'output/trimal/trimal.fas'
    log:
        'output/logs/trimal.log'
    threads:
        20
    shell:
        'bin/trimAl/source/trimal '
        '-in {input} '
        '-gt 0.6 ' ##gap threshold - fraction of sequences with gap allowed - removes all positions with gaps in more than 60% of seq.s
        '-out {output} '
        '2> {log}'

##concatenate alignments
rule FASconCAT_G_concat:
    input:
        expand('output/muscle/{gene}.fasta', gene=["ac81", "DNAP", "helicase", "lef-5",
            "lef-8", "lef-9", "p33", "pif-0", "pif-1", "pif-2", "pif-3", "pif-5"])
    output:
        'output/muscle/FcC_supermatrix.fas'
    params:
        muscle_dir = 'output/muscle'
    shell:
        'cd {params.muscle_dir} || exit 1 ; '
        'perl FASconCAT-G_v1.05.pl '
        '-s'

#################
## align genes ##
#################

##align genes separately
rule muscle_align:
    input:
        gene_fasta = 'data/core_genes/{gene}.fasta'
    output:
        alignment = 'output/muscle/{gene}.fasta'
    log:
        'output/logs/muscle_align_{gene}.log'
    threads:
        20
    shell:
        'muscle '
        '-in {input.gene_fasta} '
        '-out {output.alignment} '
        '-log {log}'

rule blast_dafv:
    input:
        dafv = 'data/DaFV.fasta'
    output:
        blastp_res = 'output/dafv_blastp/dafv_blastp.outfmt6'
    params:
        blast_db = 'bin/blast_db/nr/nr'
    threads:
        20
    log:
        'output/logs/blast_dafv.log'
    singularity:
        blast_container
    shell:
        'blastp '
        '-query {input.dafv} '
        '-db {params.blast_db} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-outfmt "6 std staxids salltitles" > {output.blastp_res} '
        '2> {log}'
