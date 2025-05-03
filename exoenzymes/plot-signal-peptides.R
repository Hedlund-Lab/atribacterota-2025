## This script generates sequence logos for the signal peptides associated with
## each of the phyla included in the analysis.
groups = c('Dictyoglomota', 'Thermotogota', 'Atribacterota', 'Verrucomicrobiota', 'Synergistota')
library('data.table')
library('magrittr')
library('ggplot2')
library('ggseqlogo')
library('cowplot')

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

signalp_sets = fread('signalp-sets.tsv.xls.gz', sep='\t', header=TRUE)
signalp_subset = signalp_sets[set != 'signalp4' & !is.na(signal) & !is.na(cs.start),][which(group == 'Other Gram -') %>% sample(1e4) %>% c(., which(group != 'Other Gram -')),]

signalp_subset[group == 'Caldatribacteriota', group := 'Atribacterota']
signalp_subset[group == 'Thermotoga', group := 'Thermotogota']
signalp_subset[group == 'Dictyoglomi', group := 'Dictyoglomota']
signalp_subset[group == 'Synergistetes', group := 'Synergistota']
signalp_subset[group == 'Verrucomicrobia', group := 'Verrucomicrobiota']

signalp_subset[set == 'intersect', set := 'SignalP 4.1 &\nSignalP 5']
signalp_subset[set == 'signalp5', set := 'SignalP 5']


dats = signalp_subset[!is.na(group),]
pdf('signal-length.pdf')
p = ggplot(dats, aes(x=set, y=cs.start)) + 
  geom_violin(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(position=position_jitter(width = 0.15, height = 0), alpha = 1, size = 0.025) +
  facet_grid(. ~ group, scales = 'free_x', space = 'free_x') +
  theme_bw() +
  xlab('Detection method + Phylum') +
  ylab('Signal peptide length')
plot(p)
dev.off()

p1 = ggplot(dats, aes(x=set, y=cs.start, fill = set)) + 
  #geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  geom_violin(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(position=position_jitter(width = 0.15, height = 0), alpha = 1, size = 0.025) +
  facet_grid(. ~ group, scales = 'free_x', space = 'free_x') +
  theme_bw()  +
  theme(strip.clip = "off") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(legend.position = 'none') +
  xlab('Detection method + Phylum') +
  ylab('Signal peptide length')


p2 = ggplot(dats, aes(x=cs.start, fill=set, y = stat(density))) + 
  geom_density(alpha = 0.5, position="identity", color = 'black') +
  facet_grid(. ~ group) +
  theme_bw() +
  theme(strip.clip = "off") +
  theme(legend.position = 'bottom') +
  xlab('Signal peptide length') +
  ylab('Density')

svg('signal-length-boxplots.svg')
    {plot_grid(p1, p2, nrow = 2) %>% plot()}
dev.off()

signalp_subset[,cs.motif.sanitized := cs.motif %>% gsub('^ +','', .) %>% gsub('-','',.)]
plots = signalp_subset[,group] %>% unique() %>% lapply(function(gi) lapply(c('SignalP 4.1 &\nSignalP 5','SignalP 5'),
    function(ii)
    {
        ggplot() +
            geom_logo(signalp_subset[set == ii & group == gi, cs.motif.sanitized], seq_type = 'aa', method = 'prob') +
            ggtitle(paste0(gi, '\n', ii)) +
            scale_x_continuous(breaks = 1:5, labels = c(-3,-2,-1,1,2)) +
            #scale_y_continuous(breaks = 0:3, labels = 0:3, limits = c(0,3.5)) +
            scale_y_continuous(breaks = 0:4/4, labels = 0:4/4, limits = c(0,1)) +
            theme_logo() +
            theme(legend.position = 'none',
                  plot.title = element_text(size=8))
            
    }    )) %>% do.call('c',.)
    
pdf('logos-prob.pdf', height = 12, width = 6)
p = plot_grid(plotlist = plots, ncol = 2)
plot(p)
dev.off()

plots = signalp_subset[,group] %>% unique() %>% lapply(function(gi) lapply(c('SignalP 4.1 &\nSignalP 5','SignalP 5'),
    function(ii)
    {
        ggplot() +
            geom_logo(signalp_subset[set == ii & group == gi, cs.motif.sanitized], seq_type = 'aa', method = 'bits') +
            ggtitle(paste0(gi, '\n', ii)) +
            scale_x_continuous(breaks = 1:5, labels = c(-3,-2,-1,1,2)) +
            scale_y_continuous(breaks = 0:3, labels = 0:3, limits = c(0,3.5)) +
            #scale_y_continuous(breaks = 0:4/4, labels = 0:4/4, limits = c(0,1)) +
            theme_logo() +
            theme(legend.position = 'none',
                  plot.title = element_text(size=8))
            
    }    )) %>% do.call('c',.)
    
pdf('logos-bits.pdf', height = 12, width = 6)
p = plot_grid(plotlist = plots, ncol = 2)
plot(p)
dev.off()

fwrite(dats[,.(cs.start, set, group)], 'signal-lengths.tsv.xls', sep='\t')

perm.len = lapply(dats[,group] %>% unique(),
function(phyl)
{
    z = dats[group == phyl & !is.na(cs.start),
            perm_test(cs.start[set == 'SignalP 5'], cs.start[set == 'SignalP 4.1 &\nSignalP 5'])
            ] %>% unclass() %>% as.list() %>% as.data.table()
    return(z[1,] %>% unlist())
}) %>% do.call('rbind',.) %>% setNames(c('len.diffs','len.pval'))
pvals = cbind(phylum = dats[,group] %>% unique(), perm.len)
colnames(pvals) = c('phylum', 'len.difference', 'len.pval')
fwrite(pvals, 'p-values-signalp5-vs-intersect.srplen.tsv.xls', sep='\t')


sink('phylum-anova-srplen.txt')
    aov(cs.start ~ group, dats[set == 'SignalP 5',]) %>% summary()
    aov(cs.start ~ group, dats[set == 'SignalP 5',]) %>% TukeyHSD()
    aov(cs.start ~ group, dats[set == 'SignalP 4.1 &\nSignalP 5',]) %>% summary()
    aov(cs.start ~ group, dats[set == 'SignalP 4.1 &\nSignalP 5',]) %>% TukeyHSD()
sink()

signalp_subset[,.(phylum = group, cs.motif.sanitized, seqID, genome, cs.start, cs.stop)] %>%
fwrite('signalp-logo-sourcedata.txt.xls', sep = '\t')