## Objective is to reproduce the dotplot in figure 3 from Katayama et al. 2020;
#  Katayama, T. et al. Isolation of a member of the candidate phylum
#  ‘Atribacteria’ reveals a unique cell membrane structure. Nat Commun 11, 6381
#  (2020).


library('data.table')
library('magrittr')
library('stringi')
library('stringr')


## SignalP 4 output: column "?" is Y for secretion and N for not.
## SignalP 5 output is tabular.
parse_signalp_outputs = function(fname)
{
    header = readLines(fname, n = 1)
    
    if (grepl('SignalP-4', header))
    {
        ## Zip out spaces and replace them with tabs, then reread the file.
        tab = fread(fname, sep='', header=FALSE)[-(1:2),V1] %>%
            stri_replace_all_regex(., ' +', '\t') %>%
            paste0(collapse='\n') %>%
            fread(sep='\t', header=FALSE)
        colnames(tab) = c('name', 'Cmax', 'pos', 'Ymax', 'pos','Smax', 'pos', 'Smean', 'D', 'prediction', 'Dmaxcut', 'Networks.used')
        
    } else if (grepl('SignalP-5', header))
    {
        tab = fread(fname, sep='\t', header=TRUE)
        colnames(tab)[1] = 'ID'
    } else {
        stop('Version information not found in file header.')
    }
    
    return(tab)
}


parse_tmhmm_outputs = function(fname)
{
    cn = c('name', 'len', 'ExpAA', 'First60', 'PredHel', 'Topology')
    type = c('character', 'numeric', 'numeric', 'numeric', 'numeric', 'character')
    tab = fread(fname, sep='\t', header=TRUE)
    colnames(tab) = cn
    
    ## Clean up and recast numeric cols.
    tab[,c(cn[-1]) := mget(cn[-1]) %>% lapply(., stri_replace_first_regex, '.*=','')]
    tab[,c(cn[type == 'numeric']) := mget(cn[type == 'numeric']) %>% lapply(., as.numeric)]
    return(tab)
}

aggregate_signalp4 = function(signalp) signalp[,sum(prediction == 'Y') / .N]
aggregate_signalp5 = function(signalp) signalp[,sum(Prediction != 'OTHER') / .N]
aggregate_tmhmm = function(tmhmm) tmhmm[,sum(Topology != 'o') / .N]

## Calc for signalp5
signalp5_dir = 'Atribacter_signalp&tmhmm/signalp5_results/'
signalp5_files = list.files(signalp5_dir, '_summary.signalp5', full.names=TRUE)
signalp5_genomes = signalp5_files %>% basename() %>% gsub('_summary.signalp5$','',.)

signalp5_dt = data.table(file = signalp5_files, genome = signalp5_genomes, version = 'signalp5')
signalp5_dt[,signalp5 := file %>% lapply(parse_signalp_outputs) %>% lapply(aggregate_signalp5) %>% unlist()]

## Calc for signalp4
signalp4_dir = 'Atribacter_signalp&tmhmm/siganlp4_results/short/'
signalp4_files = list.files(signalp4_dir, '.faa.short_out', full.names=TRUE)
signalp4_genomes = signalp4_files %>% basename() %>% gsub('.faa.short_out$','',.)

signalp4_dt = data.table(file = signalp4_files, genome = signalp4_genomes, version = 'signalp4')
signalp4_dt[,signalp4 := file %>% lapply(parse_signalp_outputs) %>% lapply(aggregate_signalp4) %>% unlist()]

## Calc for tmhmm
tmhmm_dir = 'Atribacter_signalp&tmhmm/TMHMM2/'
tmhmm_files = list.files(tmhmm_dir, '.faa.xls', full.names=TRUE)
tmhmm_genomes = tmhmm_files %>% basename() %>% gsub('.faa.xls$','',.)

tmhmm_dt = data.table(file = tmhmm_files, genome = tmhmm_genomes)
tmhmm_dt[,tmhmm := file %>% lapply(parse_tmhmm_outputs) %>% lapply(aggregate_tmhmm) %>% unlist()]


## Merge signalp
setkey(signalp4_dt, genome)
setkey(signalp5_dt, genome)

signalp_dt = signalp5_dt[signalp4_dt, .(genome, signalp5, signalp4)]
signalp_dt[,signalp_ratio := signalp5/signalp4]

## Merge in tmhmm
setkey(signalp_dt, genome)
setkey(tmhmm_dt, genome)

## Wrie to file
plotdata = signalp_dt[tmhmm_dt, .(genome, tmhmm, signalp4, signalp5, signalp_ratio)]

fwrite(plotdata, 'processed-exoenzyme-data.tsv.xls', sep='\t')

