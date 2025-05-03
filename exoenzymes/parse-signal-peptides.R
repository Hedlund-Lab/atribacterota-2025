## Parse out the signal peptides so that we can align them.
library('data.table')
library('magrittr')
library('stringi')
library('stringr')

source('io.align.R')
source('sequtils.R')

mafft = function(fname)
{
    # G-INS-i 
    cmd = paste('mafft', '--globalpair --maxiterate 1000', fname, '>', paste0(fname, '-align'))
    message(cmd)
    system(cmd)
    return(paste0(fname, '-align'))
}

hmmbuild = function(fname)
{
    
    cmd = paste('hmmbuild', paste0(fname, '.hmm'), fname)
    message(cmd)
    system(cmd)
    return(paste0(fname, '.hmm'))
}

trimal = function(fname)
{
    cmd = paste('trimal -in', fname, '-out', paste0(fname, '.trimal'), '-gappyout')
    message(cmd)
    system(cmd)
    return(paste0(fname, '.trimal'))
}

plotdata = fread('processed-exoenzyme-data.tsv.xls', sep='\t', header=TRUE)

tempura_metadata = fread('200617_TEMPURA.csv', sep=',', header=TRUE)[assembly_or_accession != '' & superkingdom != 'Archaea',]
setkey(tempura_metadata, assembly_or_accession)

gtdb_metadata = fread('op9-gtdb-metadata.tsv.xls', sep='\t', header=TRUE)
gtdb_metadata[,phylum := stri_extract_first_regex(gtdb_taxonomy, '(?<=;p__)([A-Z]|[a-z])+')]
setkey(gtdb_metadata, ncbi_genbank_assembly_accession)

plotdata[, INSDC := stri_extract_first_regex(genome, 'G[A-Z][A-Z]_[0-9]{9,9}\\.[0-9]+')]
setkey(plotdata, INSDC)

plotdata[tempura_metadata, Phylum := phylum]
plotdata[gtdb_metadata, Phylum := phylum]

plotdata[genome %in% c('OP9_SIUC_contigs', 'OP9_GBS_contigs'), Phylum := 'Caldatribacteriota']

plotdata[genome == 'OP9_SIUC_contigs', label := 'SIUC']
plotdata[genome == 'OP9_GBS_contigs', label := 'GBS']


plotdata[,group := Phylum]
plotdata[group == 'Thermotogae', group := 'Thermotogota']
plotdata[!group %in% c('Dictyoglomota', 'Thermotogota', 'Caldatribacteriota', 'Verrucomicrobia', 'Synergistetes'), group := 'Other Gram -']

plotdata[Phylum %in% c('Firmicutes', 'Actinobacteria', 'Tenericutes'), Phylum := NA]

kept = plotdata[!is.na(Phylum), genome]

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

disentangle_signalp_seqs = function(signalp4_fname, signalp5_fname)
{
    signalp4 = parse_signalp_outputs(signalp4_fname)[prediction == 'Y', name]
    signalp5fdt = parse_signalp_outputs(signalp5_fname)
    signalp5 = signalp5fdt[Prediction != 'OTHER', ID]
    csdat = signalp5fdt[Prediction != 'OTHER', get('CS Position')] %>% strsplit('\\.|: ') %>% do.call('rbind', .) %>% as.data.table()
    colnames(csdat) = c('.', 'cs', 'cs.aa', '.', '.', '.')
    csparsed = csdat[,strsplit(cs, '-') %>% do.call('rbind', .) %>% apply(2, as.numeric) %>% apply(2, as.list)] %>% cbind(signalp5, ., csdat[,cs.aa])
    colnames(csparsed) = c('signalp5id', 'cs.start', 'cs.stop', 'cs.motif')
    
    
    
    iset = intersect(signalp4, signalp5)
    sd4 = setdiff(signalp4, signalp5)
    sd5 = setdiff(signalp5, signalp4)
    out = list(data.table(seqID = iset, set = 'intersect')
        ,data.table(seqID = sd4, set = 'signalp4')
        ,data.table(seqID = sd5, set = 'signalp5'))  %>%
    do.call('rbind', .)
    
    setkey(out, seqID)
    setkey(csparsed, signalp5id)
    
    out[csparsed, cs.start := cs.start %>% unlist()] 
    out[csparsed, cs.stop := cs.stop %>% unlist()] 
    out[csparsed, cs.motif := cs.motif %>% unlist()] 
    
    return(out)
}

## Calc for signalp5
signalp5_dir = 'Atribacter_signalp&tmhmm/signalp5_results/'
signalp5_files = list.files(signalp5_dir, '_summary.signalp5', full.names=TRUE)
signalp5_genomes = signalp5_files %>% basename() %>% gsub('_summary.signalp5$','',.)
signalp5_dt = data.table(signalp5_files, genome = signalp5_genomes, version = 'signalp5')[genome %in% kept,]

## Calc for signalp4
signalp4_dir = 'Atribacter_signalp&tmhmm/siganlp4_results/short/'
signalp4_files = list.files(signalp4_dir, '.faa.short_out', full.names=TRUE)
signalp4_genomes = signalp4_files %>% basename() %>% gsub('.faa.short_out$','',.)
signalp4_dt = data.table(signalp4_files, genome = signalp4_genomes, version = 'signalp4')[genome %in% kept,]

## Merge signalp
setkey(signalp4_dt, genome)
setkey(signalp5_dt, genome)
signalp_dt = signalp5_dt[signalp4_dt, .(genome, signalp4_files, signalp5_files)]

## Pull signalp data.
signalp_sets = signalp_dt[,disentangle_signalp_seqs(signalp4_files, signalp5_files), by=genome][!is.na(seqID),]

## Get taxonomy
setkey(plotdata, genome)
setkey(signalp_sets, genome)

signalp_sets[plotdata, Phylum := Phylum]
signalp_sets[plotdata, group := group]

## Build a bunch of alignments for hmmer.
dir.create('alignments')

## Read in sequences.
faa_fnames = signalp_sets[,genome] %>% unique() %>% paste0('faa/', ., '.faa')
faa = lapply(faa_fnames, read.align, truncate_at = -1) %>% do.call('c', .)

signalp_sets[,signal := substr(faa[seqID] %>% unlist(), 1, cs.start)]
signalp_sets[,full_sequence := faa[seqID] %>% unlist()]

fwrite(signalp_sets, 'signalp-sets.tsv.xls.gz', sep='\t')
#groups = signalp_sets[Phylum == group, group] %>% unique()

dir.create('signal-alignments')
dir.create('full-seqs')
#groups = c('Dictyoglomota', 'Thermotogota', 'Caldatribacteriota', 'Verrucomicrobia', 'Synergistetes')
# groups = c('Verrucomicrobia', 'Synergistetes')
# for(gi in groups)
# {
    # fname_signalp5 = paste0('signal-alignments/', gi, '-signalp5.faa')
    # fname_intersect = paste0('signal-alignments/', gi, '-intersect.faa')
    # gi_signalp_sets = signalp_sets[group == gi,]
    # gi_signalp_sets[set == 'signalp5', as.list(signal)] %>% setNames(gi_signalp_sets[set == 'signalp5',seqID]) %>% write.align(fname_signalp5)
    # gi_signalp_sets[set == 'intersect', as.list(signal)] %>% setNames(gi_signalp_sets[set == 'intersect',seqID]) %>% write.align(fname_intersect)
    
    
    # fname_mafft_signalp5 = fname_signalp5 %>% mafft()
    # fname_mafft_signalp5 %>% hmmbuild()
    # fname_mafft_signalp5 %>% trimal () %>% hmmbuild()
    
    # fname_mafft_intersect = fname_intersect %>% mafft()
    # fname_mafft_insersect %>% hmmbuild()
    # fname_mafft_insersect %>% trimal() %>% hmmbuild()
# }

groups = c('Dictyoglomota', 'Thermotogota', 'Caldatribacteriota', 'Verrucomicrobia', 'Synergistetes')
# groups = c('Verrucomicrobia', 'Synergistetes')
for(gi in groups)
{
    fname_signalp5 = paste0('full-seqs/', gi, '-signalp5.faa')
    fname_intersect = paste0('full-seqs/', gi, '-intersect.faa')
    gi_signalp_sets = signalp_sets[group == gi,]
    gi_signalp_sets[set == 'signalp5', as.list(full_sequence)] %>% setNames(gi_signalp_sets[set == 'signalp5',seqID]) %>% write.align(fname_signalp5)
    gi_signalp_sets[set == 'intersect', as.list(full_sequence)] %>% setNames(gi_signalp_sets[set == 'intersect',seqID]) %>% write.align(fname_intersect)
    
    
    # fname_mafft_signalp5 = fname_signalp5 %>% mafft()
    # fname_mafft_signalp5 %>% hmmbuild()
    # fname_mafft_signalp5 %>% trimal () %>% hmmbuild()
    
    # fname_mafft_intersect = fname_intersect %>% mafft()
    # fname_mafft_insersect %>% hmmbuild()
    # fname_mafft_insersect %>% trimal() %>% hmmbuild()
}



# alignment = read.align('signal-alignments/Caldatribacteriota-intersect.faa-align') %>% unlist()

# ggplot() +
    # geom_logo(alignment, seq_type = 'aa', method = 'bits') +
    # #scale_y_continuous(breaks = 0:4/4, labels = 0:4/4, limits = c(0,1)) +
    # theme_logo() +
    # theme(legend.position = 'none',
          # plot.title = element_text(size=8))

