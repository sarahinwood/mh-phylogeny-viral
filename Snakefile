#!/usr/bin/env python3

blast_container = 'docker://ncbi/blast:2.12.0'
raxmlng_container = 'docker://evolbioinfo/raxml-ng:v1.0.2'

#########
# RULES #
#########

rule target:
    input:
        #'output/blastp_viral/blastp_viral.outfmt6',
        'output/dafv_blastp/dafv_blastp.outfmt6',
        'output/raxml/tree_inf/01.raxml.bestTree',
        'output/raxml/tree_inf/02.raxml.bestTree',
        'output/raxml/raxml_all/DCMut.raxml.bestTree'

rule raxml_all:
    input:
        msa = 'output/raxml/check/check.raxml.reduced.phy'
    output:
        'output/raxml/raxml_all/DCMut.raxml.bestTree'
    params:
        wd = 'output/raxml/raxml_all',
        prefix = 'output/raxml/raxml_all/DCMut'
    singularity:
        raxmlng_container
    log:
        'output/logs/raxml_all.log'
    threads:
        20
    shell:
        'raxml-ng '
        '--msa {input.msa} '
        '--model DCMut ' ##check this is best model, change number bootstrap reps
        '--prefix {params.prefix} '
        '--threads {threads} '
        '--seed 2 '
        '2> {log}'

rule raxml_search1:
    input:
        msa = 'output/raxml/check/check.raxml.reduced.phy'
    output:
        'output/raxml/tree_inf/02.raxml.bestTree'
    params:
        wd = 'output/raxml/tree_inf',
        prefix = 'output/raxml/tree_inf/02'
    singularity:
        raxmlng_container
    log:
        'output/logs/raxml_tree_inf_search.log'
    threads:
        20
    shell:
        'raxml-ng '
        '--search1 '
        '--msa {input.msa} '
        '--model LG '
        '--prefix {params.prefix} '
        '--threads {threads} '
        '--seed 2 '
        '2> {log}'

rule raxml_tree_inf:
    input:
        msa = 'output/raxml/check/check.raxml.reduced.phy'
    output:
        'output/raxml/tree_inf/01.raxml.bestTree'
    params:
        wd = 'output/raxml/tree_inf',
        prefix = 'output/raxml/tree_inf/01'
    singularity:
        raxmlng_container
    log:
        'output/logs/raxml_tree_inf.log'
    threads:
        20
    shell:
        'mkdir {params.wd} || exit 1 ;'
        'raxml-ng '
        '--msa {input.msa} '
        '--model LG '
        '--prefix {params.prefix} '
        '--threads {threads} '
        '--seed 2 '
        '--tree pars{{25}},rand{{25}} ' ##25 parsimony trees and 25 random trees
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
rule trimal:
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

rule blastp_viral:
    input:
        prodigal_proteins = 'data/Mh_prodigal/protein_translations.faa',
        blast_db = 'output/blastdb/viral_db/viral_db.phr'
    output:
        blastp_res = 'output/blastp_viral/blastp_viral.outfmt6'
    params:
        viral_db = 'output/blastdb/viral_db/viral_db'
    threads:
        20
    log:
        'output/logs/blastp_viral.log'
    shell:
        'blastp '
        '-query {input.prodigal_proteins} '
        '-db {params.viral_db} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-outfmt "6 std salltitles" > {output.blastp_res} '
        '2>{log}'

rule make_viral_blast_db:
    input:
        'output/viral_proteins.fasta'
    output:
        blast_db = 'output/blastdb/viral_db/viral_db.phr'
    params:
        db_name = 'viral_db',
        db_dir = 'output/blastdb/viral_db/viral_db'
    log:
        'output/logs/make_viral_blast_db.log'
    shell:
        'makeblastdb '
        '-in {input} '
        '-dbtype prot '
        '-title {params.db_name} '
        '-out {params.db_dir} '
        '-parse_seqids '
        '2> {log}'

rule cat_viral_files:
    output:
        'output/viral_proteins.fasta'
    shell:
        'cat data/ncbi_proteins/*.fasta > {output}'