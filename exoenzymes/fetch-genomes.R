## Fetch outgroup genomes
library('data.table')
library('magrittr')

source('io.align.R')

gtdb_ref_data_path = '/mnt/d/Ubuntu/science/reference/gtdb/release202/'
taxonomy_path = paste0(gtdb_ref_data_path, '/taxonomy/', 'gtdb_taxonomy.tsv')

taxonomy = fread(taxonomy_path, sep='\t', header=FALSE) %>% setkey(V1)

metadata_path = paste0(gtdb_ref_data_path, '/metadata/', 'bac120_metadata_r202.tsv.gz')
metadata = fread(metadata_path, sep='\t', header=TRUE)

op9s = taxonomy[ grepl('^d__Bacteria;p__Caldatribacteriota;', V2)
               | grepl('^d__Bacteria;p__Dictyoglomota;', V2)
               | grepl('^d__Bacteria;p__Thermotogota;', V2)
               
               
               , V1]
op9md = metadata[accession %in% op9s,]

op9hq = op9md[checkm_completeness >= 90 & checkm_contamination <= 5,]

filepaths = fread(paste0(gtdb_ref_data_path, '/fastani/genome_paths.tsv'), sep=' ', header=FALSE) %>% setkey(V1)

op9hq[,unique(accession)] %>% lapply(
    function(x)
    {
        if (!dir.exists(paste0('contigs/',x)))
        {
            filename = paste0(x,'_genomic.fna.gz') %>% gsub('^.._','',.)
            filepath = filepaths[filename, V2]
            cmd = paste0('gunzip -ckq ', gtdb_ref_data_path, 'fastani/', filepath, filename,' > contigs/',x,'.fna')
            message(paste0('Copying ', x, '...'))
            message(cmd)
            system(cmd)
        }
    }) %>% invisible()


fwrite(op9hq, 'op9-gtdb-metadata.tsv.xls', sep='\t')


prodigal = function(contig, gff, faa = NULL, cds = NULL)
{
    flag_gff = paste('-o', gff)
    flag_faa = ifelse(is.null(faa), '', paste('-a', faa))
    flag_cds = ifelse(is.null(cds), '', paste('-d', cds))
    flag_contig = paste('-i', contig)
    cmd = paste('prodigal', flag_gff, flag_contig, flag_faa, flag_cds, '-q')
    message(cmd)
    system(cmd)
}

contigs = list.files('contigs/', '.fna$')

dir.create('faa')
dir.create('cds')
dir.create('gff')

for (cf in contigs)
{
    contig = paste0('contigs/', cf)
    gn = gsub('\\.fna$', '', cf)
    
    faa = paste0('./faa/', gn, '.faa')
    cds = paste0('./cds/', gn, '.fca')
    gff = paste0('./gff/', gn, '.gff')
    if (file.exists(contig)) prodigal(contig, gff, faa, cds)
}
