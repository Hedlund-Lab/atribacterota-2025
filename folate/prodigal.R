library('data.table')
library('magrittr')

taxa = list.files('fna', '\\.fna$') %>% gsub('\\.fna$', '', .)

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

for (ti in taxa)
{
    dir.create(paste0('./faa/', ti), recursive = TRUE, showWarnings = FALSE)
    dir.create(paste0('./cds/', ti), recursive = TRUE, showWarnings = FALSE)
    dir.create(paste0('./gff/', ti), recursive = TRUE, showWarnings = FALSE)    
    contig = paste0('fna/', ti, '.fna')
    faa = paste0('./faa/', ti,'/', ti , '.faa')
    cds = paste0('./cds/', ti,'/', ti , '.fca')
    gff = paste0('./gff/', ti,'/', ti , '.gff')
    if (file.exists(contig) &
          (
                !file.exists(faa)
           |    !file.exists(cds)
           |    !file.exists(gff)
          )
    ) prodigal(contig, gff, faa, cds)
}

for (ti in taxa)
{
    cmd = paste('cat', paste0(paste0('./cds/', ti, '/*.fca')), '>', paste0('cds/', ti, '.fca'))
    message(cmd)
    system(cmd)
    
    cmd = paste('cat', paste0(paste0('./faa/', ti, '/*.faa')), '>', paste0('faa/', ti, '.faa'))
    message(cmd)
    system(cmd)
}