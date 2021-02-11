README for annotations/

Records of table origins:

##############
### RefSeq ###
##############

For refseq gene set from hg19, "hg19_RefSeq_UCSC.gtf":

http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=570932845_m2xM6kxkc5VU6aNFrBguBIMVbcOC&clade=mammal&org=Human&db=hg19&hgta_group=genes&hgta_track=refGene&hgta_table=0&hgta_regionType=genome&position=chr21%3A33031597-33041570&hgta_outputType=primaryTable&hgta_outFileName=

see "hg19_RefSeq_UCSC_request.png" for options. Last accessed 2016/12/08.


#############################
### Transcription Factors ###
#############################

For a list of all human transcription factors, see the following most-inclusive set of TFs:

A census of human transcription factors: function, expression and evolution
Juan M Vaquerizas, Sarah K Kummerfeld, Sarah A Teichmann, Nicholas M Luscombe
Nature Reviews Genetics 10, 252-263 (April 2009) | doi:10.1038/nrg2538

- Supplementary Table 2: http://www.nature.com/nrg/journal/v10/n4/extref/nrg2538-s3.txt downloaded on 2017/01/04. See “annotations/TFtable.txt”.

Opened in Excel and copied a column 2 (Ensembl ID) to a new file for simple manipulation in R. See annotations/TF_gene_ids.csv


Also, on The Feb 8, downloaded Human_TF_MasterList_v1_02.xlsx from http://humantfs.ccbr.utoronto.ca/download/Human_TF_MasterList_v1_02.xlsx. This is a newer census of human transcription factors with citation:
Lambert SA, Jolma A, Campitelli LF, Das PK, Yin Y, Albu M, Chen X, Taipale J, Hughes TR, Weirauch MT.(2018) The Human Transcription Factors. Cell. 172(4):650-665. doi: 10.1016/j.cell.2018.01.029. Review.

Opened in Excel, copied columns [3:end, A:H], renamed colnames, and saved as \t-delimited txt file “annotations/2018_HTFreview.txt” for manipulation in R.



#############
### HTSeq ###
#############

hg19_forHTSeq.gtf used for HTSeq-0.6.1 count analysis. From …




################
### gene_ids ###
################

hg19gene_idToGeneSymbol.tsv used in downstream analysis from …



#######################
### gene_id lengths ###
#######################

hg19_lengths_per_gene_id.tsv calculated from hg19_forHTseq.gtf using GenomicFeatures package v1.26.3. Full commands:
txdb <- makeTxDbFromGFF("~/Dropbox (RajLab)/Projects/cellid/annotations/hg19_forHTSeq.gtf", format="gtf")
lengthsPergeneid <- sum(width(IRanges::reduce(exonsBy(txdb, by = "gene")))) # a named int(), in bases
lengthtbl<-as_tibble(list(
  gene_id = names(lengthsPergeneid),
  length = lengthsPergeneid
))
write.table(lengthtbl, file = "~/Dropbox (RajLab)/Projects/cellid/annotations/hg19_lengths_per_gene_id.tsv", sep = "\t", quote = F, row.names = F)



