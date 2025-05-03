library('data.table')
library('magrittr')

fbcogs = fread('folate-biosynthesis-cog.txt', sep='\t', header=TRUE) %>% setkey(cog)
fbkos  = fread('folate-biosynthesis-ko.txt', sep='\t', header=FALSE) %>% setkey(V1)


kos  = fread('ko.tsv.xls', sep='\t', header=TRUE)
cogs = fread('cogs.parsed.txt', sep='\t', header=TRUE)

dats = rbind(cogs[COG %in% fbcogs[,cog],.(gene, annot=COG)], kos[ko %in% fbkos[,V1], .(gene=seq, annot=ko)])[order(gene),]
dats[annot %in% fbkos[,V1], desc := fbkos[annot,V2]]
dats[annot %in% fbcogs[,cog], desc := fbcogs[annot,desc]]

fwrite(dats, 'cogs-kos-folate-combined.tsv.xls', sep='\t')