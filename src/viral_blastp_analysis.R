library(data.table)
library(dplyr)
library(rtracklayer)

viral_blast <- fread("output/blastp_viral/blastp_viral.outfmt6")
gff <- readGFF("data/Mh_prodigal/gene_predictions.gff")
gff$gene_no <- tstrsplit(gff$ID, "_", keep=c(2))

##viral_contig_blast res
setnames(viral_blast, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"),
         new=c("peptide_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(viral_blast, peptide_id, evalue, -bit_score)
##make species column, remove [] around name
viral_blast$species <- gsub(".+\\[(.+)\\]", "\\1", viral_blast$annotation)
viral_blast$annotation <- tstrsplit(viral_blast$annotation, " \\[", keep=c(1))

##table of species name and counts - all results
all_species <- count(viral_blast, species, sort=TRUE)

##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
min_evalues <- viral_blast[,.SD[which.min(evalue)], by=peptide_id]
##table of species name and counts - best results
best_species <- count(min_evalues, species, sort=TRUE)

##want to see which peptides had hits to multiple viruses
##to decide which genes to use for phylogeny
peptide_by_virus <- viral_blast[,c(1, 14)]
peptide_by_virus_un <- unique(peptide_by_virus)
peptides_species_counts <- count(peptide_by_virus_un, peptide_id, sort=TRUE)
sum(peptides_species_counts$n>3)

##peptides & counts by best hit
min_evalue_annot <- min_evalues[,c(1,3,11,13,14)]
full_table <- merge(peptides_species_counts, min_evalue_annot, by="peptide_id", all.x=TRUE)
no_bro <- subset(full_table, !grepl("bro", full_table$annotation))
no_bro <- subset(no_bro, !grepl("baculovirus repeat ORF", no_bro$annotation))
sum(no_bro$n>3)

fwrite(no_bro, "output/blastp_viral/annots_peptide_species_counts.csv")

