## Parse hmmer file and output to matching TSVs.
## Input format: hmmer domtblout and/or hmmer tblout from hmmsearch, and a
## string specifying the mode (domtblout or tblout)
## Output format: per-gene TSV showing whole-alignment and per-domain hits.

args = commandArgs(trailingOnly = TRUE)
infile = args[1]
mode = args[2]

parse_and_output_hmmer = function(fname, mode, sort.column = NULL)
{
    txt = data.table::fread(fname, header = FALSE, data.table=FALSE, sep='', skip = 3)[,1]
    txt = txt[-1*(length(txt):(length(txt)-9))]
    txtsplit = vapply(txt, strsplit, split = ' +', list(1))
    if(mode == 'tblout')
    {
        txtmerge = do.call('rbind',unname(lapply(txtsplit, function(x) c(x[1:18],paste(x[19:length(x)], collapse= ' ')))))
        colnames(txtmerge) = c('target name', 'target accession', 'query name', 'query accession', 'full sequence E-value',
                               'full sequence score', 'full sequence bias', 'best domain E-value', 'best domain score',
                               'best domain bias', 'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'description of target')
    } else if (mode == 'domtblout') {
        txtmerge = do.call('rbind',unname(lapply(txtsplit, function(x) c(x[1:22],paste(x[23:length(x)], collapse= ' ')))))
        colnames(txtmerge) = c('target name', 'target accession', 'tlen', 'query name', 'query accession', 'qlen',
                               'full sequence E-value', 'full sequence score', 'full sequence bias', 'this domain #',
                               'this domain of', 'this domain c-Evalue', 'this domain i-Evalue', 'this domain score',
                               'this domain bias', 'hmm from', 'hmm coord to', 'ali from', 'ali coord to', 'env from',
                               'env coord to', 'acc', 'description of target')
    }
    
    if (!is.null(sort.column))
    {
        txtmerge = txtmerge[order(txtmerge[,sort.column]),]
    }
    
    data.table::fwrite(txtmerge, paste0(fname, '.ps.tsv'), row.names=FALSE, quote=FALSE, sep='\t')
    return()
}

parse_and_output_hmmer(infile, mode, 1)