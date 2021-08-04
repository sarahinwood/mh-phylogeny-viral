library(data.table)
library(dplyr)

dafv_blast <- fread("output/dafv_blastp/dafv_blastp.outfmt6")

##viral_contig_blast res
setnames(dafv_blast, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14"),
         new=c("peptide_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "taxid", "annotation"))
##remove hits to self
nonself <- subset(dafv_blast, !(grepl("Drosophila-associated", annotation)))

##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(nonself, peptide_id, evalue, -bit_score)
min_evalues <- nonself[,.SD[which.min(evalue)], by=peptide_id]
