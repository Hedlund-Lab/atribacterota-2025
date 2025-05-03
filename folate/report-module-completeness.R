library('magrittr')
library('data.table')
library('stringi')
library('openxlsx')

#gtdb_ref_data_path = '/mnt/d/Ubuntu/science/reference/gtdb/release207/'
#taxonomy_path = 'gtdb-taxonomy-atribacterota.tsv'

#taxtxt = fread('atribacterota-genomes.txt', sep='\t', header=TRUE) %>% setkey(genome)

functionName = 'vitamin-biosynthesis'

keggReferenceDir = '/mnt/d/Ubuntu/science/reference/ko-2024-07-08/'

keggProfilesDir = paste0(keggReferenceDir, '/profiles')
moduleDefinitionsFile = 'module_data.tsv'
moduleDefinitionsFilepath = paste0(keggReferenceDir, moduleDefinitionsFile)

moduleFile = 'vitamin-biosynthesis-module.txt'
koFile = 'vitamin-biosynthesis-ko.txt'

moduleDefinitions = fread(moduleDefinitionsFilepath, sep='\t', header=FALSE) %>% setkey(V1)

module = fread(moduleFile, sep='\t', header=FALSE)
kos = fread(koFile, sep='\t', header=FALSE)
module[,searchTerm := paste0('md:', V1)]

setkey(module, searchTerm)

selectModuleDefinitions = moduleDefinitions[module, .(V1, V2, V3)]

submodules = selectModuleDefinitions[, .(M = strsplit(V3, "(\\((?:[^()]++|(?1))*\\))(*SKIP)(*F)|( |\\+)", perl=TRUE)[[1]]),by=.(md=V1)]

submodules[,M := stringi::stri_replace_all_fixed(M, ' \\-\\-','')]
submodules[,M := stringi::stri_replace_all_regex(M, '\\-K[0-9][0-9][0-9][0-9][0-9]','')]

submodules[,M := stringi::stri_replace_all_fixed(M, ' ',' & ')]
submodules[,M := stringi::stri_replace_all_fixed(M, '+',' & ')]
submodules[,M := stringi::stri_replace_all_fixed(M, ',',' | ')]
submodules[,weight := 1]
submodules = submodules[M!='--',]

moduleKO = selectModuleDefinitions[,.(KO = V3 %>% stri_extract_all_regex(., 'K[0-9]{5,5}') %>% unlist() %>% unique()), by=.(Module=V1)] %>% setkey(KO)

annotations = fread('ko.tsv.xls', sep='\t', header=TRUE) 
annotations[,genome := stri_replace_first_regex(seq, '\\|.*','')]
annotations %>% setkey(ko, genome)

genomeLookup = rbind(annotations[,.(inko = genome %>% unique(), acc = genome %>% unique())],
                     annotations[,.(inko = genome %>% unique(), acc = genome %>% unique() %>% paste0('GB',.))],
                     annotations[,.(inko = genome %>% unique(), acc = genome %>% unique() %>% gsub('GCA','GCF',.) %>% paste0('RS',.))]) %>% setkey(acc)

intaxa = taxtxt[genomeLookup, .(taxonomy = gtdb_taxonomy, genome = inko)][!is.na(taxonomy),]
taxonomy = intaxa[,strsplit(taxonomy, '; ')[[1]], by=genome][,strsplit(V1,'__') %>% do.call('rbind',.) %>% as.data.table(), by=genome] %>% dcast(genome ~ V1, value.var='V2')
taxonomy[, d := gsub('^.$','', d)]
taxonomy[, p := gsub('^.$','', p)]
taxonomy[, c := gsub('^.$','', c)]
taxonomy[, o := gsub('^.$','', o)]
taxonomy[, f := gsub('^.$','', f)]
taxonomy[, g := gsub('^.$','', g)]
taxonomy[, s := gsub('^.$','', s)]

phyla = taxonomy[order(d, p, c, o, f, g, s),.(genome, phylum = p)]

orderedphyla = rbind(phyla[phylum == 'Atribacterota',],phyla[phylum != 'Atribacterota',])[phylum != '',]

genomes = orderedphyla[,genome]
genomes = genomes[!genomes %in% c('GCA_011370655.1', 'GCF_008630935.1')] %>% gsub('^GCA','_GCA',.)

moduleKOGrid = moduleKO[,.(m = rep(Module, length(genomes)), k = rep(KO, length(genomes)), g = lapply(genomes, rep, .N) %>% unlist())] %>% setkey(k, g)
KOGrid = moduleKO[,.(k = rep(KO, length(genomes)), g = lapply(genomes, rep, .N) %>% unlist())] %>% unique()
KOGrid[,g := factor(g, levels = genomes)]

KOGrid = KOGrid %>% setkey(k, g)


koGrid = annotations[KOGrid, .(genome, seq, KO = k)][order(KO),] %>% unique()
koGrid[, n := (!is.na(seq)) * 1]
koMatrix = koGrid %>% dcast(genome ~ KO, fun.aggregate = sum, value.var = 'n') %>% na.omit()

submoduleCompleteness = submodules[,M] %>% lapply(function(x) with(koMatrix, eval(parse(text=x)) > 0) %>% as.numeric()) %>% do.call('rbind',.) %>% as.data.table() %>% setNames(koMatrix[,genome]) %>% cbind(.,submodules[,.(md,M)])
submoduleMelt = submoduleCompleteness %>% melt()
moduleCompleteness = submoduleMelt[,sum(value)/length(value),by=.(genome = variable, module = md)]

moduleCompletenessMatrix = moduleCompleteness %>% dcast(genome ~ module, value.var = 'V1') %>% setkey(genome)


moduleGrid = annotations[moduleKOGrid, .(genome, seq, Module = m, KO = k, ko)][order(Module, KO),] %>% unique()
moduleGrid[, n := (!is.na(seq)) * 1]
moduleGrid[, code := paste0(Module, '-', KO)]
moduleKOMatrix = moduleGrid %>% dcast(genome ~ code, fun.aggregate = sum, value.var = 'n') %>% na.omit() %>% setkey(genome)

nomodule = annotations[!ko %in% moduleKO[,KO] & ko %in% kos[,V1],][,table(ko) %>% as.data.table(), by=genome][,.(genome, ko = paste0('NoModule-',ko), N)] %>% dcast(genome ~ ko, value.var = 'N')
nomodule %>% setkey(genome)
unorderedMatrix = nomodule[moduleKOMatrix,][moduleCompletenessMatrix,]
unorderedMatrix[is.na(unorderedMatrix)] = 0
unorderedMatrix[,NoModule := '-']
unorderedMatrix %>% setkey(genome)

orderedMatrix = unorderedMatrix[genomes,-1][,mget(unorderedMatrix[,-1] %>% colnames() %>% sort())] %>% cbind(unorderedMatrix[genomes,1],.)

wb = createWorkbook()
addWorksheet(wb, "KEGG Modules")

## Decorate table.
## rule applies to all each cell in range
zeroStyle = createStyle(fontColour = "#CCCCCC", bgFill = "#FFFFFF",valign='center', halign='center')
oneStyle = createStyle(bgFill = "#CCCCCC", fontColour = "#000000",valign='center', halign='center')
twoStyle = createStyle(bgFill = "#AAAAAA", fontColour = "#000000",valign='center', halign='center')
moduleStyle = createStyle(border=c("top", "bottom", "left", "right"), textDecoration='bold')
speciesStyle = createStyle(fgFill = "#FFFFFF", fontColour = "#000000", textDecoration='bold', valign='center', halign='right')
phylumStyle = createStyle(fgFill = "#FFFFFF", fontColour = "#000000", textDecoration='bold', valign='top', halign='left')

innerData = orderedMatrix[,-1]
moduleCols = which(grepl('M[0-9]{5,5}$',colnames(innerData)) | colnames(innerData) == 'NoModule')
nmoduleCols = which(colnames(innerData) == 'NoModule')
koCols = setdiff(1:ncol(innerData), moduleCols)

# Subheader
shStyle = createStyle(textRotation = 90, valign='bottom',halign='center')
subheadings = matrix(gsub('.*[-:]','',gsub('-$','',colnames(innerData))), nrow=1)
writeData(wb, "KEGG Modules", subheadings, startCol=3, startRow=2, colNames = FALSE)
addStyle(wb, "KEGG Modules", shStyle, cols=(1:ncol(subheadings))+2, rows=2, gridExpand=TRUE)
setRowHeights(wb, "KEGG Modules", rows = 2, heights = 57)
setRowHeights(wb, "KEGG Modules", rows = 1, heights = 51)

# Header
membership = cumsum(1:ncol(innerData) %in% moduleCols)
headingModules = subheadings[1,moduleCols[membership]]

mData = selectModuleDefinitions
mData[,m := gsub('.*\\:','',V1)]
setkey(mData, m)
mData[,desc := paste0(m, ' : ', V2)]

headingLabels = matrix(mData[headingModules, desc], nrow=1)
headingLabels[is.na(headingLabels)] = 'Unassigned'
writeData(wb, "KEGG Modules", headingLabels, startCol=3, startRow=1, colNames = FALSE)

mergeUs = data.table(group=membership, ii=1:ncol(innerData))[,.(MIN=min(ii), MAX=max(ii)), by=group]
for(cc in 1:nrow(mergeUs))
{
    mergeCells(wb, "KEGG Modules", mergeUs[cc,c(MIN,MAX)+2], 1)
}

## Sidebar.
op = taxonomy[,.(taxon = paste(d, p, o, f, g, s, sep=';')), by=genome]
op %>% setkey(genome)
orderedphyla %>% setkey(genome)
writeData(wb, "KEGG Modules", op[orderedMatrix[,genome],taxon] %>% matrix(ncol=1), startCol=1, startRow=3, colNames = FALSE)
membership2 = orderedphyla[orderedMatrix[,genome],cumsum(!duplicated(phylum))]
mergeUs2 = data.table(group=membership2, ii=1:nrow(innerData))[,.(MIN=min(ii), MAX=max(ii)), by=group]
for(cc in 1:nrow(mergeUs2))
{
    #mergeCells(wb, "KEGG Modules", rows = mergeUs2[cc,c(MIN,MAX)+2], cols = 1)
    #addStyle(wb, "KEGG Modules", phylumStyle, rows = mergeUs2[cc,(MIN:MAX)+2], cols = 1)
}

rows = matrix(orderedMatrix[,genome], ncol=1)
writeData(wb, "KEGG Modules", rows, startCol=2, startRow=3, colNames = FALSE)
setColWidths(wb, "KEGG Modules", cols=1, widths = 21)
setColWidths(wb, "KEGG Modules", cols=2, widths = 21)


topStyle = createStyle(border="top", borderStyle='medium')
bottomStyle = createStyle(border="bottom", borderStyle='medium')

addStyle(wb, "KEGG Modules", topStyle, cols=(1:(nrow(innerData)+2)), rows=mergeUs2[,MIN+2], gridExpand=TRUE)
addStyle(wb, "KEGG Modules", bottomStyle, cols=(1:(nrow(innerData)+2)), rows=mergeUs2[,MAX+2], gridExpand=TRUE)

## Draw the body of the table.
writeData(wb, "KEGG Modules", innerData, startCol=3, startRow=3, colNames = FALSE)
conditionalFormatting(wb, "KEGG Modules", cols=(koCols)+2, rows=(1:nrow(innerData))+2, rule="=0", style = zeroStyle)
conditionalFormatting(wb, "KEGG Modules", cols=(koCols)+2, rows=(1:nrow(innerData))+2, rule="=1", style = oneStyle)
conditionalFormatting(wb, "KEGG Modules", cols=(koCols)+2, rows=(1:nrow(innerData))+2, rule=">1", style = twoStyle)
conditionalFormatting(wb, "KEGG Modules", cols=(moduleCols)+2, rows=(1:nrow(innerData))+2, rule="<0.5", style = zeroStyle)
conditionalFormatting(wb, "KEGG Modules", cols=(moduleCols)+2, rows=(1:nrow(innerData))+2, rule=">=0.5", style = oneStyle)
conditionalFormatting(wb, "KEGG Modules", cols=(moduleCols)+2, rows=(1:nrow(innerData))+2, rule=">=0.8", style = twoStyle)
addStyle(wb, "KEGG Modules", moduleStyle, cols=(moduleCols)+2, rows=(1:nrow(innerData))+2, gridExpand=TRUE)
setColWidths(wb, "KEGG Modules", cols=(1:ncol(innerData))+2, widths = ifelse(1:ncol(innerData) %in% koCols, 1.5, 3.57),
    hidden = c(colSums(innerData[,1:(nmoduleCols-1)]), 1, colSums(innerData[,(nmoduleCols+2):ncol(innerData)]))  == 0)


freezePane(wb, "KEGG Modules", firstActiveRow = 3, firstActiveCol = 3)

writeData(wb, "KEGG Modules", matrix("Module Definitions", nrow=1), startCol=1, startRow=(nrow(innerData))+4, colNames=FALSE)
writeData(wb, "KEGG Modules", mData, startCol=1, startRow=(nrow(innerData))+5, colNames=FALSE)

unlink("module-completeness-report.xlsx")
saveWorkbook(wb, file = "module-completeness-report.xlsx", overwrite = TRUE)