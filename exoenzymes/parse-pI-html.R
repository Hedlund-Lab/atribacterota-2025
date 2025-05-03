## Parse out pI html files.
library('data.table')
library('magrittr')
library('stringi')

source('io.align.R')

hydrophobicity = fread('hydrophobicity.tsv.xls', sep='\t', header=TRUE)
hydrophobicityScale = hydrophobicity[,setNames(kd.Hydrophobicity, abbr)]

## get pI and hydrophobicity
pIdt = list.files('signal-alignments/','-pI.html') %>% lapply(
function(f)
{
    html = paste0('signal-alignments/', f) %>% fread(sep='', header=FALSE)
    IPC = html[,V1] %>% stri_extract_all_regex(., 'IPC peptide \\(pI=[^(]+\\)') %>% unlist() %>% na.omit() %>% gsub('IPC peptide (pI=', '', ., fixed=TRUE) %>% gsub(').*', '', .) %>% as.numeric()
    GENE = html[,V1] %>% stri_extract_all_regex(., 'Input sequence:<br><p>&gt;.+$') %>% unlist() %>% na.omit() %>% gsub('Input sequence:<br><p>&gt;', '', ., fixed=TRUE)
    
    data.table(f, IPC, GENE)
}) %>% do.call('rbind', .) %>% setkey(GENE)


hydrophobicitydt = list.files('signal-alignments/','.faa$') %>% lapply(
function(f)
{
    faa = paste0('signal-alignments/', f) %>% read.align()
    GENE = names(faa)
    HYD = faa %>% unlist() %>% strsplit('') %>% lapply(function(x) hydrophobicityScale[x]) %>% lapply(mean) %>% unlist()
    
    data.table(f, HYD, GENE) %>% na.omit()
}) %>% do.call('rbind', .) %>% setkey(GENE)

dats = hydrophobicitydt[pIdt, .(phylum = gsub('-.*','',f), set = gsub('.*-','',f) %>% gsub('\\.faa$','',.), pI = IPC, hydrophobicity = HYD, GENE)] %>% na.omit()

fwrite(dats, 'hydrophobicity-and-pI.tsv.xls', sep='\t')