library('data.table')
library('magrittr')
library('ggplot2')
library('ggrepel')
library('showtext')
library('stringi')
library('lme4')

logit = function(p) log(p / (1 - p));

pal = c('Atribacterota' = '#F12029', 'Dictyoglomota' = '#F4961D', 'Thermotogota' = '#36B249','Other Gram -' = '#A3A3A3', 'Gram +' = '#333333')
plotdata = fread('processed-exoenzyme-data.tsv.xls', sep='\t', header=TRUE)

tempura_metadata = fread('200617_TEMPURA.csv', sep=',', header=TRUE)[assembly_or_accession != '' & superkingdom != 'Archaea',]
setkey(tempura_metadata, assembly_or_accession)

gtdb_metadata = fread('op9-gtdb-metadata.tsv.xls', sep='\t', header=TRUE)
gtdb_metadata[,phylum := stri_extract_first_regex(gtdb_taxonomy, '(?<=;p__)([A-Z]|[a-z])+')]
setkey(gtdb_metadata, ncbi_genbank_assembly_accession)

plotdata[, INSDC := stri_extract_first_regex(genome, 'G[A-Z][A-Z]_[0-9]{9,9}\\.[0-9]+')]
setkey(plotdata, INSDC)

plotdata[tempura_metadata, Phylum := phylum]
plotdata[tempura_metadata, GC := Genome_GC]
plotdata[tempura_metadata, Topt := Topt_ave]
plotdata[tempura_metadata, Genome_size := Genome_size]
plotdata[gtdb_metadata, Phylum := phylum]

plotdata[genome %in% c('OP9_SIUC_contigs', 'OP9_GBS_contigs'), Phylum := 'Atribacterota']

plotdata[genome == 'OP9_SIUC_contigs', label := 'SIUC']
plotdata[genome == 'OP9_GBS_contigs', label := 'GBS']

plotdata[Phylum == 'Caldatribacteriota', group := 'Atribacterota']
plotdata[Phylum == 'Caldatribacteriota', Phylum := 'Atribacterota']
plotdata[,group := Phylum]
plotdata[group == 'Thermotogae', group := 'Thermotogota']
plotdata[group == 'Thermotoga', group := 'Thermotogota']
plotdata[group == 'Dictyoglomi', group := 'Dictyoglomota']
plotdata[group == 'Synergistetes', group := 'Synergistota']
plotdata[group == 'Verrucomicrobia', group := 'Verrucomicrobiota']
plotdata[!group %in% c('Dictyoglomota', 'Thermotogota', 'Atribacterota'), group := 'Other Gram -']

plotdata[Phylum %in% c('Firmicutes', 'Actinobacteria', 'Tenericutes'), group := 'Gram +']
#plotdata[Phylum %in% c('Firmicutes', 'Actinobacteria', 'Tenericutes'), Phylum := NA]

plotdata[,group := factor(group, levels = names(pal))]

## Report stats
sink('exoenzyme-statistics.txt')
    summary(lm(log2(signalp_ratio) ~ group, data = plotdata)) %>% print() ## Log transform signalp ratios
    summary(lm(logit(tmhmm) ~ group, data = plotdata)) %>% print() ## Logit transform percent of genome with TM domains
sink()

p = ggplot(data=plotdata[!is.na(Phylum),]) +
    geom_point(aes(x = tmhmm*100, y = signalp_ratio, color=group)) +
    stat_ellipse(aes(x = tmhmm*100, y = signalp_ratio, color=group, fill = group), geom = 'polygon', alpha = 0.15, linetype = 'dashed') +
    geom_text_repel(aes(x = tmhmm*100, y = signalp_ratio, label = label), family = 'Arial') +
    scale_color_manual(breaks = names(pal), values = pal) +
    scale_fill_manual(breaks = names(pal), values = pal) +
    theme_bw() +
    xlab('Proteins encoding transmembrane helices\n(% of total genes)') +
    ylab('Proteins encoding Sec signal peptides\n(ratio between results from Signalp 5.0 / Signalp 4.1)') +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(family = 'Arial'),
      axis.title.x = element_text(family = 'Arial', face = 'bold'),
      axis.text.y = element_text(family = 'Arial'),
      axis.title.y = element_text(family = 'Arial', face = 'bold'),
      legend.text = element_text(family = 'Arial'),
      legend.title = element_text(family = 'Arial', face = 'bold'),
      )

svg('exoenzyme-dotplot.svg', width = 7, height = 6)
    plot(p)
dev.off()