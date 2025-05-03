library('data.table')
library('dplyr')
library('stringi')

trim_leading_whitespace = function(x) x %>% stri_replace_first_regex('^ +', '')
trim_trailing_whitespace = function(x) x %>% stri_replace_last_regex(' +$', '')
trim_outside_whitespace = function(x) x %>% trim_leading_whitespace() %>% trim_trailing_whitespace()
reduce_whitespace = function(x) x %>% stri_replace_all_regex(' +',' ')


args = commandArgs(trailingOnly = TRUE)
infile = args[1]

txt = fread(infile, sep='\t', header=FALSE, colClasses='character')[!is.na(V1),]

txt[,isNewQuery := stri_detect_fixed(V1, 'Query=')]
txt[,isEndQuery := nchar(V1) == 0]
txt[,querySxn := cumsum(isNewQuery)]
txt[, collapsed := ifelse(any(isNewQuery) & any(isEndQuery), paste0(V1[which(isNewQuery):which(isEndQuery)], collapse=''), ''), by=querySxn]


filtered = txt[substr(V1, 1, 1) == '>' | isNewQuery,]
rm(txt)
invisible(gc())

unparsedAnnots = filtered[,.(collapsed[isNewQuery][1], V1[!isNewQuery]), by=querySxn][!is.na(V2),]

parsedAnnots = unparsedAnnots[,.(gene = stri_replace_first_fixed(V1, 'Query=','') %>% trim_outside_whitespace() %>% gsub(' .*','',.), V2)][,.(COG = stri_extract_all_regex(V2, 'COG[0-9]{4,4}') %>% unlist()) , by=gene] %>% unique()
boundAnnots = parsedAnnots[!duplicated(gene),]

outfile = gsub('\\.[^.]+$', '.parsed.txt',infile)
fwrite(boundAnnots, outfile, sep='\t', quote=FALSE)