library('data.table')
library('magrittr')
library('ggplot2')
library('cowplot')
#library('kSamples')

## permutation test
perm_test = function(x, y, n_perm = 1e3)
{
  obsdiff = mean(x) - mean(y)
  
  ## Calc differences after resampling
  diffs = data.table(perm = 1:n_perm)[,
  {
    resample = c(x, y) %>% sample()
    .(diff = mean(resample[1:length(x)]) - mean(resample[-(1:length(x))]))
  }, by = perm][,diff]
  
  pval = mean(abs(diffs) >= abs(obsdiff))
  
  return(list(obsdiff, pval))
}

dats = fread('hydrophobicity-and-pI.tsv.xls', sep='\t', header=TRUE)

dats[set == 'intersect', set := 'SignalP 4.1 &\nSignalP 5']
dats[set == 'signalp5', set := 'SignalP 5']

dats[phylum == 'Caldatribacteriota', phylum := 'Atribacterota']
dats[phylum == 'Thermotoga', phylum := 'Thermotogota']
dats[phylum == 'Dictyoglomi', phylum := 'Dictyoglomota']
dats[phylum == 'Synergistetes', phylum := 'Synergistota']
dats[phylum == 'Verrucomicrobia', phylum := 'Verrucomicrobiota']

plots = dats[,phylum] %>% unique() %>% lapply(
function(gi)
    {
        ggplot() +
            geom_point(data = dats[phylum == gi,], aes(x = hydrophobicity, y = pI, color = set), alpha = 0.2) +
            ggtitle(gi) +
            theme_bw() +
            #scale_y_continuous(breaks = 0:4/4, labels = 0:4/4, limits = c(0,1)) +
            theme(legend.position = 'none',
                  plot.title = element_text(size=8))
    })

plots2 = c(plots, list(get_legend(plots[[length(plots)]] + theme(legend.position = 'right'))))

pdf('hydrophobicity-and-pI.pdf', height = 12, width = 6)
    p = plot_grid(plotlist = plots2, ncol = 2)
    plot(p)
dev.off()

#pdf('isoelectric-point.pdf')
p1 = ggplot(dats, aes(x=set, y=pI, fill = set)) + 
  #geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  #geom_jitter(position=position_jitter(width = 0.3, height = 0), alpha = 0.2, aes(color = set)) +
  geom_violin(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(position=position_jitter(width = 0.15, height = 0), alpha = 1, size = 0.025) +
  facet_grid(. ~ phylum, scales = 'free_x', space = 'free_x') +
  theme_bw()  +
  theme(strip.clip = "off") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(legend.position = 'none') +
  xlab('Detection method + Phylum') +
  ylab('Isoelectric point')
#plot(p1)
#dev.off()

#pdf('hydrophobicity.pdf')
p2 = ggplot(dats, aes(x=set, y=hydrophobicity, fill = set)) + 
  #geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  #geom_jitter(position=position_jitter(width = 0.3, height = 0), alpha = 0.2, aes(color = set)) +
  geom_violin(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(position=position_jitter(width = 0.15, height = 0), alpha = 1, size = 0.025) +
  facet_grid(. ~ phylum, scales = 'free_x', space = 'free_x') +
  theme_bw()  +
  theme(strip.clip = "off") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(legend.position = 'none') +
  xlab('Detection method + Phylum') +
  ylab('Hydrophobicity')
#plot(p2)

#dev.off()


#pdf('hydrophobicity-hist.pdf')
p3 = ggplot(dats, aes(x=hydrophobicity, fill=set, y = stat(density))) + 
  geom_histogram(alpha = 0.5, position="identity") +
  facet_grid(. ~ phylum) +
  theme_bw() +
  theme(strip.clip = "off") +
  theme(legend.position = 'bottom') +
  xlab('Hydrophobicity') +
  ylab('Density')
#plot(p3)

#dev.off()

#pdf('pI-hist.pdf')
p4 = ggplot(dats, aes(x=pI, fill=set, y = stat(density))) + 
  geom_histogram(alpha = 0.5, position="identity") +
  facet_grid(. ~ phylum) +
  theme_bw() +
  theme(strip.clip = "off") +
  theme(legend.position = 'bottom') +
  xlab('pI') +
  ylab('Density')
#plot(p4)

#dev.off()

svg('pI.svg')
    {plot_grid(p1, p4, nrow = 2) %>% plot()}
dev.off()

svg('hydrophobicity.svg')
    {plot_grid(p2, p3, nrow = 2) %>% plot()}
dev.off()

svg('pI-hydrophobicity-boxplots.svg')
    {plot_grid(p1, p2, nrow = 2) %>% plot()}
dev.off()

perm.pI = lapply(dats[,phylum] %>% unique(),
function(phyl)
{
    z = with(dats[phylum == phyl,], 
        perm_test(pI[set == 'SignalP 5'], pI[set == 'SignalP 4.1 &\nSignalP 5'])) %>% unclass() %>% as.list() %>% as.data.table()
    return(z[1,] %>% unlist())
}) %>% do.call('rbind',.) %>% setNames(c('pI.diffs','pI.pval'))

perm.hydrophobicity = lapply(dats[,phylum] %>% unique(),
function(phyl)
{
    z = with(dats[phylum == phyl,], 
        perm_test(hydrophobicity[set == 'SignalP 5'], hydrophobicity[set == 'SignalP 4.1 &\nSignalP 5'])) %>% unclass() %>% as.list() %>% as.data.table()
    return(z[1,] %>% unlist())
}) %>% do.call('rbind',.)  %>% setNames(c('hydrophobicity.diffs','hydrophobicity.pval'))

pvals = cbind(phylum = dats[,phylum] %>% unique(), perm.hydrophobicity, perm.pI)
colnames(pvals) = c('phylum', 'hydrobicity.difference', 'hydrohobicity.pval', 'pI.difference', 'pI.pval')
fwrite(pvals, 'p-values-signalp5-vs-intersect.tsv.xls', sep='\t')

sink('phylum-anova.txt')
    aov(hydrophobicity ~ phylum, dats[set == 'SignalP 5',]) %>% summary()
    aov(hydrophobicity ~ phylum, dats[set == 'SignalP 5',]) %>% TukeyHSD()
    aov(hydrophobicity ~ phylum, dats[set == 'SignalP 4.1 &\nSignalP 5',]) %>% summary()
    aov(hydrophobicity ~ phylum, dats[set == 'SignalP 4.1 &\nSignalP 5',]) %>% TukeyHSD()
    
    aov(pI ~ phylum, dats[set == 'SignalP 5',]) %>% summary()
    aov(pI ~ phylum, dats[set == 'SignalP 5',]) %>% TukeyHSD()
    aov(pI ~ phylum, dats[set == 'SignalP 4.1 &\nSignalP 5',]) %>% summary()
    aov(pI ~ phylum, dats[set == 'SignalP 4.1 &\nSignalP 5',]) %>% TukeyHSD()
sink()
