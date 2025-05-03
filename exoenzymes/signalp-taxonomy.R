## Process taxonomy information for SignalP entries in taxa of interest.
source('io.align.R')

library('data.table')
library('magrittr')
library('stringi')

gtdb_ref_data_path = '/mnt/d/Ubuntu/science/reference/gtdb/release202/'
taxonomy_path = paste0(gtdb_ref_data_path, '/taxonomy/', 'gtdb_taxonomy.tsv')

taxonomy = fread(taxonomy_path, sep='\t', header=FALSE) %>% setkey(V1)

taxonomy[,V2 := gsub('Verrucomicrobiota','Verrucomicrobia',V2)]
taxonomy[,V2 := gsub('Synergistota','Synergistetes',V2)]


taxoncounts = 
c('Caldatribacteriota','Thermotogota','Dictyoglomota', 'Verrucomicrobia', 'Synergistetes') %>% lapply(
function(x)
{
    cols = c('qseqid','sseqid','pident','length', 'mismatch','gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
    signalp5 = paste0('blout/', x, '-signalp5.blout.tsv') %>% fread(sep='\t') %>% setNames(cols)
    signalp5[,sgenome := gsub('\\|.*','', sseqid)]
    signalp5_len = paste0('signal-alignments/', x, '-signalp5.faa') %>% read.align()  %>% length() %>% setNames('signalp5')


    intrsect = paste0('blout/', x, '-intersect.blout.tsv') %>% fread(sep='\t') %>% setNames(cols)
    intrsect[,sgenome := gsub('\\|.*','', sseqid)]
    intersect_len = paste0('signal-alignments/', x, '-intersect.faa') %>% read.align() %>% length() %>% setNames('intersect')


    blout = rbind(intrsect[,.(qseqid, sgenome, bitscore, set = 'intersect')], signalp5[,.(qseqid, sgenome, bitscore, set = 'signalp5')])
    setkey(blout, sgenome)

    blout[taxonomy, phylum := stri_extract_first_regex(V2, '(?<=;p__)([A-Z]|[a-z])+')]

    nblout = blout[phylum != x,]
    maxblout = nblout[sgenome %in% nblout[,.(sgenome[bitscore == max(bitscore)][1]),by=.(qseqid)][,V1],]

    bloutcounts = maxblout[,table(phylum) %>% as.data.table(), by = set]
    bloutcounts[,taxon := x]
    bloutcounts[set == 'intersect', P := N/intersect_len]
    bloutcounts[set == 'signalp5', P := N/signalp5_len]
    
    return(bloutcounts)
}) %>% do.call('rbind',.)

taxonall = rbind(taxoncounts, taxoncounts[,.(N = NA, phylum = 'Unassigned', P = 1 - sum(P)), by = .(set, taxon)])

fwrite(taxonall, 'signal-peptide-taxonomy-proportions.tsv.xls', sep='\t')