#######################################
# The following script includes analysis for the
# gavaging gnotobiotic mice with in vitro communities
# that were derived from the small intestine and stool 
# of a single healthy, human subject
# 
# Rebecca Culver
#######################################

# configure directories, load libraries and base functions
rm(list=ls())
source(paste0(here::here(), "/0-config.R"))
ps<-readRDS(paste0(clean_data_dir,'merge_phyloseq.RDS'))
genus.order<-readRDS(paste0(fig_dir, 'official_colors.RDS'))

#######################################
# Fig. 2A: In vitro communities used in mouse experiments
# Includes assessment of their similarity to the human stool and small intestine sample as well
#######################################

# We will merge all replicates of StoolCom and SICom to get the most accurate picture of how
# representative the in vitro community is to the human inoculum
ps2plot<- ps %>%
  subset_samples(Project %in% 'MouseExp2' | SampleType %in% c('Overnight','original_2522_capsule',"Human_Stool")) %>% 
  merge_samples("Community")

taxa2keep<-psmelt(ps2plot %>% transform_sample_counts(function(x) x/sum(x))) %>%
  filter(Abundance > 0.001) 

ps2plot <- prune_taxa(taxa2keep$OTU, ps2plot) %>% transform_sample_counts(function(x) x/sum(x))

df2plot <- psmelt(ps2plot) %>% 
  mutate(Sample = factor(Sample, levels=c('Human_Intestinal_Microbiota','IC','Human_Stool','SC'))) %>%
  mutate(Family = factor(ifelse(is.na(Family), ifelse(!is.na(Phylum), paste0("(Unclassified ", Phylum, ")"), "(Unclassified phylum)"), Family))) %>%
  left_join(genus.order, by="Family") %>%
  mutate(Family = factor(Family, levels=genus.order$Family)) %>%
  group_by(Sample, Family, colors) %>%
  dplyr::summarise(Abundance=sum(Abundance)) %>% ungroup() 

(a<-ggplot(df2plot , aes(x=Sample, y=Abundance, fill=Family)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=df2plot$colors) +
    theme_minimal() +
    theme(panel.margin.x=unit(0, "lines")) +
    theme(axis.text.y= element_text(color='black'))+ #, axis.text.x=element_blank(), axis.ticks.x=element_blank())+   
    ylab('Rel. ab.')   +           
    xlab('') )                                                                     

ggsave(paste0(fig_dir, 'subpanels/Fig_2A_barplots.pdf'), width = 7, height = 4)

### Similarity to human inoculum
df.similarity<-psmelt(ps2plot) %>% 
  filter(Abundance > 0.001) %>%
  # Re-label NA families
  mutate(Family = factor(ifelse(is.na(Family), ifelse(!is.na(Phylum), paste0("(Unclassified ", Phylum, ")"), "(Unclassified phylum)"), Family))) 

# ASV count
df.similarity %>% group_by(Sample) %>% dplyr::summarise(n())

# More dataframes by sample
df.stool <- df.similarity %>% filter(Sample %in% 'Human_Stool')
df.stoolcom <- df.similarity %>% filter(Sample %in% 'SC')
df.intestine <- df.similarity %>% filter(Sample %in% 'Human_Intestinal_Microbiota')
df.sicom <- df.similarity %>% filter(Sample %in% 'IC')

# Family abundance analysis
a<-unique(df.stool$Family)
length(a)
b<-unique(df.stoolcom$Family)
length(b)
intersect(a,b)
df.stool %>% filter(Family %in% df.stoolcom$Family) %>% dplyr::summarise(Abund=sum(Abundance))

a<-unique(df.intestine$Family)
length(a)
b<-unique(df.sicom$Family)
length(b)
intersect(a,b)
df.intestine %>% filter(Family %in% df.sicom$Family) %>% dplyr::summarise(Abund=sum(Abundance))

# Correlation analysis - ASV level
# The text includes the non-bonferroni corrected p-value, as it was the statistic for that one correlation
# rather than a p value associated as significant across a collection of correlations

ps.corr <- ps2plot %>%
  transform_sample_counts(fun=filterfun1) %>% # turns all 0s into 0.001
  transform_sample_counts( function(x)  log10(x))

asvs2correlate <- otu_table(ps.corr)
rownames(asvs2correlate) <- rownames(sample_data(ps.corr))

res <- rcorr(t(asvs2correlate), type="pearson") # compute correlations
corr.df.tmp <- flattenCorrMatrix(res$r, res$P) 
p.vals <- p.adjust(corr.df.tmp$p ,method="bonferroni") # adjust pvalues
r.vals <- corr.df.tmp$cor
corr.df <- cbind(corr.df.tmp %>% select(-p), p.vals)
corr.df

# Correlation analysis - Family level
# Because certain families appear as NA, we must rename the tax table dataframe
ps.corr.fam <- ps2plot

taxa.renamed <- data.frame(tax_table(ps.corr.fam)) %>%
  mutate(Family = factor(ifelse(is.na(Family), ifelse(!is.na(Phylum), paste0("(Unclassified ", Phylum, ")"), "(Unclassified phylum)"), Family)))
tax_table(ps.corr.fam) <- as.matrix(taxa.renamed)

ps.corr.fam <- ps.corr.fam %>%
  tax_glom("Family") %>%
  transform_sample_counts(fun=filterfun1) %>% # turns all 0s into 0.001
  transform_sample_counts( function(x)  log10(x))

asvs2correlate <- otu_table(ps.corr.fam)
rownames(asvs2correlate) <- rownames(sample_data(ps.corr.fam))

res <- rcorr(t(asvs2correlate), type="pearson") # compute correlations
corr.df.tmp <- flattenCorrMatrix(res$r, res$P) 
p.vals <- p.adjust(corr.df.tmp$p ,method="bonferroni") # adjust pvalues
r.vals <- corr.df.tmp$cor
corr.df <- cbind(corr.df.tmp %>% select(-p), p.vals)
corr.df
#######################################
# Fig. 2B: Mouse experiment overview
# Schematic created in illustrator
#######################################

#######################################
# Fig. 2C: Relative abundance of taxa in each mouse
#######################################

# Mouse samples dataframe
mouse.ps <- ps %>% subset_samples(Project %in% 'mouse_exp2' & Location %in% c('Stool','Small Intestine')) %>%
  merge_samples("mouse_location")
taxa2keep <- mouse.ps %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% filter(Abundance > 0.001)
mouse.df <- prune_taxa(taxa2keep$OTU, mouse.ps) %>% transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% 
  mutate(Mouse = str_split(Sample, '_', simplify=TRUE)[,1],
         Location = str_split(Sample, '_', simplify=TRUE)[,2],
         Diet = str_split(Sample, '_', simplify=TRUE)[,3],
         Community = str_split(Sample, '_', simplify=TRUE)[,4],
         MouseType = str_split(Sample, '_', simplify=TRUE)[,5],
         mouse_subj = paste0(Mouse,'_',MouseType,'_',Community),
         com_loc = paste0(Community,'_',Location))

# Now add in colors and sum by family abundance
mouse.df.barplots <-  mouse.df %>%
  mutate(Diet = ifelse(Diet %in% 'Normal','Normal','MD')) %>%
  mutate(Diet = factor(Diet, levels=c('Normal','MD'))) %>%
  mutate(Community = ifelse(Community %in% 'SC','StoolCom','DualCom')) %>%
  mutate(Family = factor(ifelse(is.na(Family), ifelse(!is.na(Phylum), paste0("(Unclassified ", Phylum, ")"), "(Unclassified phylum)"), Family))) %>%
  left_join(genus.order, by="Family") %>%
  mutate(Family = factor(Family, levels=genus.order$Family)) %>%
  group_by(Sample, Location, Diet, Community, Family, colors) %>%
  dplyr::summarise(Abundance=sum(Abundance))

(c<-ggplot(mouse.df.barplots , aes(x=Sample, y=Abundance, fill=Family)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=mouse.df.barplots$colors) +
    theme_minimal() +
    theme(panel.margin.x=unit(0, "lines")) +
    theme(axis.text.y= element_text(color='black'), axis.text.x=element_blank(), axis.ticks.x=element_blank())+                                                                                                                                            
    ylab('Rel. ab.')   +                                                                                                                                                                                                                                                                                         
    facet_wrap(~Location+Community+Diet, scale='free_x',nrow=2))     

ggsave(paste0(fig_dir, 'subpanels/Fig_2C_barplots.pdf'), width = 5, height = 3)

#######################################
# Supplementary Table 1: Stool and small intestine samples are highly correlated within a cage 
#######################################

mouse.ps <- ps %>% subset_samples(Project %in% 'mouse_exp2' & Location %in% c('Stool','Small Intestine')) %>%
  merge_samples("mouse_sampleType")
taxa2keep <- mouse.ps %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% filter(Abundance > 0.001)

ps.corr <- prune_taxa(taxa2keep$OTU, mouse.ps) %>% 
  transform_sample_counts(function(x) x/sum(x))  %>%
  transform_sample_counts(fun=filterfun1) %>% # turns all 0s into 0.001
  transform_sample_counts( function(x)  log10(x))

asvs2correlate <- otu_table(ps.corr)
rownames(asvs2correlate) <- rownames(sample_data(ps.corr))

res <- rcorr(t(asvs2correlate), type="pearson") # compute correlations
corr.df.tmp <- flattenCorrMatrix(res$r, res$P) 
p.vals.adj <- p.adjust(corr.df.tmp$p ,method="bonferroni") # adjust pvalues
r.vals <- corr.df.tmp$cor
corr.df <- cbind(corr.df.tmp %>% select(-p), p.vals.adj) %>%
  dplyr::rename(Sample1=row, Sample2=column) %>%
  mutate(Mouse1 = str_split(Sample1, '_', simplify=TRUE)[,1]) %>%
  mutate(Mouse2 = str_split(Sample2, '_', simplify=TRUE)[,1]) %>%
  mutate(Location = str_split(Sample1, '_', simplify=TRUE)[,2]) %>%
  mutate(Location2 = str_split(Sample2, '_', simplify=TRUE)[,2]) %>%
  mutate(Diet = str_split(Sample1, '_', simplify=TRUE)[,3]) %>%
  mutate(Community = str_split(Sample1, '_', simplify=TRUE)[,4]) %>%
  mutate(Diet2 = str_split(Sample2, '_', simplify=TRUE)[,3]) %>%
  mutate(Community2 = str_split(Sample2, '_', simplify=TRUE)[,4]) %>%
  mutate(compare_me = ifelse(Location %in% 'Stool' & Location2 %in% 'Stool',1,
                             ifelse(Mouse1==Mouse2 & !Location %in% 'Stool' & !Location2 %in% 'Stool',1,0))) %>%
  filter(Diet==Diet2, Community==Community2, compare_me==1) %>%
  select(-Community2,  -Diet2, -Sample1, -Sample2, -compare_me) %>%
  mutate(Community = ifelse(Community %in% 'SC', 'StoolCom','DualCom')) %>%
  mutate(Diet = ifelse(Diet %in% 'noFiber','MD','SD')) %>%
  arrange(Location, Community, Diet)
write.csv(corr.df, paste0(fig_dir, 'Supplementary Table 1.csv'))


#######################################
# Fig. 2D: Analysis of taxa that colonized the mouse intestinal tract
# NORMAL DIET ONLY
#######################################

# Emergent ASVs analysis:

# Get ASVs/families in the in vitro communities

overnight.cultures <- ps %>% subset_samples(Project %in% 'MouseExp2' | SampleType %in% 'Overnight') %>% 
  merge_samples("Community")
taxa2keep <- overnight.cultures %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% filter(Abundance > 0.001)
dual.com <- prune_taxa(taxa2keep$OTU, overnight.cultures) %>% transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% filter(Abundance > 0.001)
unique(dual.com$Family)
stool.com <- prune_taxa(taxa2keep$OTU, overnight.cultures) %>% transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% filter(Sample %in% 'SC', Abundance > 0.001)
unique(stool.com$Family)

#Mouse samples dataframe
mouse.ps <- ps %>% subset_samples(Project %in% 'mouse_exp2' & !SampleType %in% 'Cecum') %>%
  merge_samples("mouse_location")
taxa2keep <- mouse.ps %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% filter(Abundance > 0.001)
mouse.df <- prune_taxa(taxa2keep$OTU, mouse.ps) %>% transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% filter(Abundance > 0.001) %>%
  mutate(Mouse = str_split(Sample, '_', simplify=TRUE)[,1],
         Location = str_split(Sample, '_', simplify=TRUE)[,2],
         Diet = str_split(Sample, '_', simplify=TRUE)[,3],
         Community = str_split(Sample, '_', simplify=TRUE)[,4],
         MouseType = str_split(Sample, '_', simplify=TRUE)[,5],
         mouse_subj = paste0(Mouse,'_',MouseType,'_',Community),
         com_loc = paste0(Community,'_',Location)) %>%
  filter(Diet %in% 'Normal') # only looking at normal diet for now

dim(mouse.df)
#[1] 933  41

# StoolCom Mice
stoolcom.si <- mouse.df %>% filter(com_loc=='SC_Small Intestine') %>%
  filter(OTU %in% stool.com$OTU) 
unique(stoolcom.si$Family)

print('emergent ASVS in the small intestine of StoolCom-colonized mice')
stoolcom.si <- mouse.df %>% filter(com_loc=='SC_Small Intestine') %>%
  filter(!OTU %in% stool.com$OTU) 
stoolcom.si %>% group_by(Sample) %>% dplyr::summarise(Abund=sum(Abundance)) %>% dplyr::summarise(mean(Abund))
# Barplot
df2plot<-stoolcom.si %>% left_join(genus.order, by=c('Family')) %>%
  group_by(mouse_subj, Family, colors) %>%
  dplyr::summarise(Abundance=sum(Abundance)) %>%
  ungroup() %>% mutate(Family = factor(Family, levels=genus.order$Family)) 
cols<-df2plot %>% data.frame() %>% select(Family, colors) %>% unique() %>% mutate(Family = factor(Family, levels=genus.order$Family)) %>% arrange(Family)
(stoolcom.si.emergent.barplot <- ggplot(df2plot, aes(x=mouse_subj, y=Abundance, fill=Family)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=cols$colors) +
    theme_minimal() +
    ylim(0,1) +
    theme(panel.margin.x=unit(0, "lines")) +
    theme(axis.text.y= element_text(color='black'), axis.text.x=element_blank(), axis.ticks.x=element_blank(),,legend.position='none')+                                                                                                                                            
    ylab('Rel. ab.')  )


print('emergent ASVS in the stool of StoolCom-colonized mice')
stoolcom.stool <- mouse.df %>% filter(com_loc=='SC_Stool') %>%
  filter(OTU %in% stool.com$OTU) 
unique(stoolcom.stool$Family)
stoolcom.stool <- mouse.df %>% filter(com_loc=='SC_Stool') %>%
  filter(!OTU %in% stool.com$OTU) 
stoolcom.stool %>% group_by(Sample) %>% dplyr::summarise(Abund=sum(Abundance)) %>% dplyr::summarise(mean(Abund))
# Barplot
df2plot<-stoolcom.stool %>% left_join(genus.order, by=c('Family')) %>%
  group_by(mouse_subj, Family, colors) %>%
  dplyr::summarise(Abundance=sum(Abundance)) %>%
  ungroup() %>% mutate(Family = factor(Family, levels=genus.order$Family)) 
cols<-df2plot %>% data.frame() %>% select(Family, colors) %>% unique() %>% mutate(Family = factor(Family, levels=genus.order$Family)) %>% arrange(Family)
(stoolcom.stool.emergent.barplot <- ggplot(df2plot, aes(x=mouse_subj, y=Abundance, fill=Family)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=cols$colors) +
  theme_minimal() +
  ylim(0,1)+
  theme(panel.margin.x=unit(0, "lines")) +
  theme(axis.text.y= element_text(color='black'), axis.text.x=element_blank(), axis.ticks.x=element_blank(),,legend.position='none')+                                                                                                                                            
  ylab('Rel. ab.')  )


# DualCom Mice
dualcom.si <- mouse.df %>% filter(com_loc=='Stool&SI_Small Intestine') %>%
  filter(OTU %in% dual.com$OTU) 
unique(dualcom.si$Family)
dualcom.si <- mouse.df %>% filter(com_loc=='Stool&SI_Small Intestine') %>%
  filter(!OTU %in% dual.com$OTU) 
dualcom.si %>% group_by(Sample) %>% dplyr::summarise(Abund=sum(Abundance)) %>% dplyr::summarise(mean(Abund))
# Barplot
df2plot<-dualcom.si %>% left_join(genus.order, by=c('Family')) %>%
  group_by(mouse_subj, Family, colors) %>%
  dplyr::summarise(Abundance=sum(Abundance)) %>%
  ungroup() %>% mutate(Family = factor(Family, levels=genus.order$Family)) 
cols<-df2plot %>% data.frame() %>% select(Family, colors) %>% unique() %>% mutate(Family = factor(Family, levels=genus.order$Family)) %>% arrange(Family)
(dualcom.si.emergent.barplot <- ggplot(df2plot, aes(x=mouse_subj, y=Abundance, fill=Family)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=cols$colors) +
    theme_minimal() +
    ylim(0,1)+
    theme(panel.margin.x=unit(0, "lines")) +
    theme(axis.text.y= element_text(color='black'), axis.text.x=element_blank(), axis.ticks.x=element_blank(),legend.position='none')+                                                                                                                                            
    ylab('Rel. ab.')  )


dualcom.stool <- mouse.df %>% filter(com_loc=='Stool&SI_Stool') %>%
  filter(OTU %in% dual.com$OTU) 
unique(dualcom.stool$Family)
dualcom.stool <- mouse.df %>% filter(com_loc=='Stool&SI_Stool') %>%
  filter(!OTU %in% dual.com$OTU) 
dualcom.stool %>% group_by(Sample) %>% dplyr::summarise(Abund=sum(Abundance)) %>% dplyr::summarise(mean(Abund))
# Barplot
df2plot<-dualcom.stool %>% left_join(genus.order, by=c('Family')) %>%
  group_by(mouse_subj, Family, colors) %>%
  dplyr::summarise(Abundance=sum(Abundance)) %>%
  ungroup() %>% mutate(Family = factor(Family, levels=genus.order$Family)) 
cols<-df2plot %>% data.frame() %>% select(Family, colors) %>% unique() %>% mutate(Family = factor(Family, levels=genus.order$Family)) %>% arrange(Family)
(dualcom.stool.emergent.barplot <- ggplot(df2plot, aes(x=mouse_subj, y=Abundance, fill=Family)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=cols$colors) +
    theme_minimal() +
    ylim(0,1)+
    theme(panel.margin.x=unit(0, "lines")) +
    theme(axis.text.y= element_text(color='black'), axis.text.x=element_blank(), axis.ticks.x=element_blank(),legend.position='none')+                                                                                                                                            
    ylab('Rel. ab.')  )


plot_grid(stoolcom.si.emergent.barplot,dualcom.si.emergent.barplot, 
          stoolcom.stool.emergent.barplot, dualcom.stool.emergent.barplot, nrow = 2)

ggsave(paste0(fig_dir, 'subpanels/Fig_2D_emergentBarplots.pdf'),height=4, width=4)

#######################################
# Fig. 2E: Emergent ASV correlation plots
#######################################

# Correlation plots
# Get physeq object and merge all mice
mouse.on.ps <- ps %>% subset_samples(Project %in% 'mouse_exp2' & !SampleType %in% c('Cecum','Control'))
sample_data(mouse.on.ps) <- data.frame(mouse.on.ps@sam_data) %>%
  mutate(rep_merged= paste0(Location,'_',Diet,'_',Community))
mouse.on.ps <- mouse.on.ps %>% merge_samples("rep_merged") 

# Get the overnight communities
overnight.cultures <- ps %>% subset_samples(Project %in% 'MouseExp2' | SampleType %in% 'Overnight') %>% 
  merge_samples("Community")
taxa2keep <- overnight.cultures %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% filter(Abundance > 0.001)
dual.com <- prune_taxa(taxa2keep$OTU, overnight.cultures) %>% merge_samples("finalDaySamples") %>% transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% filter(finalDaySamples==0) %>% mutate(Sample = "DualCom")
stool.com <- prune_taxa(taxa2keep$OTU, overnight.cultures) %>% transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% filter(Sample %in% 'SC')

mouse.on.df <- mouse.on.ps %>% transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% select(-rep_merged) %>%
  rbind(dual.com, stool.com) %>%
  mutate(Family = factor(ifelse(is.na(Family), ifelse(!is.na(Phylum), paste0("(Unclassified ", Phylum, ")"), "(Unclassified phylum)"), Family))) %>%
  select(OTU, Sample, Abundance, Phylum, Family, Genus, Species) %>%
  ##group_by(Sample, Family) %>%
  #dplyr::summarise(Abundance = sum(Abundance)) %>% 
  filter(Abundance > 0) %>%
  mutate(Abundance = log10(Abundance+0.001)) %>%
  select(OTU,Sample, Family, Abundance) %>%
  pivot_wider(values_from = Abundance, names_from = Sample, values_fill=-3) %>% data.frame() %>%
  left_join(genus.order, by="Family") %>%
  mutate(Family = factor(Family, levels=genus.order$Family))

colnames(mouse.on.df)

a<-ggplot(mouse.on.df,
          aes(x=SC, y=Small.Intestine_Normal_SC)) +
  geom_point(color=mouse.on.df$colors) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  coord_fixed(1/1) +
  scale_y_continuous(breaks=c(-1,-2,-3))+
  scale_x_continuous(breaks=c(-1,-2,-3))+
  #stat_cor(method='pearson') +
  xlab('StoolCom')+
  ylab('Small intestine of StoolCom-colonized mice')

b<-ggplot(mouse.on.df,
          aes(x=SC, y=Stool_Normal_SC)) +
  geom_point(color=mouse.on.df$colors) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  coord_fixed(1/1) +
  scale_y_continuous(breaks=c(-1,-2,-3))+
  scale_x_continuous(breaks=c(-1,-2,-3))+
  #stat_cor(method='pearson') +
  xlab('StoolCom')+
  ylab('Stool of StoolCom-colonized mice')


c<-ggplot(mouse.on.df,
       aes(x=DualCom, y=Small.Intestine_Normal_Stool.SI)) +
  geom_point(color=mouse.on.df$colors) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  coord_fixed(1/1) +
  scale_y_continuous(breaks=c(-1,-2,-3))+
  scale_x_continuous(breaks=c(-1,-2,-3))+
  #stat_cor(method='pearson') +
  xlab('DualCom')+
  ylab('Small intestine of DualCom-colonized mice')

d<-ggplot(mouse.on.df,
       aes(x=DualCom, y=Stool_Normal_Stool.SI)) +
  geom_point(color=mouse.on.df$colors) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  coord_fixed(1/1) +
  scale_y_continuous(breaks=c(-1,-2,-3))+
  scale_x_continuous(breaks=c(-1,-2,-3))+
  #stat_cor(method='pearson') +
  xlab('DualCom')+
  ylab('Stool of DualCom-colonized mice')

draft<-plot_grid(plot_grid(a,b),
                 plot_grid(c,d), nrow=2, rel_heights = c(1,1))
ggsave(paste0(fig_dir, 'subpanels/Fig_2E_scatterPlots.pdf'), width = 5, height = 5)

#######################################
# Fig. 2F Shannon diversity
#######################################

# Shannon diversity
# Create data frame where each sample has corresponding Shannon div
ps.shannon<- mouse.ps <- ps %>% subset_samples(Project %in% 'mouse_exp2' & Location %in% c('Stool','Small Intestine')) %>%
  merge_samples("mouse_location") %>%
  prune_taxa(taxa_sums(.) > 50, .) # because we are inserting information about raw reads, we do want to apply some sort of filter
# Since we filter for >5000 reads, we'll say taxa sums must be at least 50 (0.1%)
# Though without the filter, the data does not change

Sample<-rownames(ps.shannon@sam_data)

sample_data(ps.shannon) <- data.frame(ps.shannon@sam_data) %>% 
  cbind(Sample) %>%
  mutate(Mouse = str_split(Sample, '_', simplify=TRUE)[,1],
         Location = str_split(Sample, '_', simplify=TRUE)[,2],
         Diet = str_split(Sample, '_', simplify=TRUE)[,3],
         Community = str_split(Sample, '_', simplify=TRUE)[,4],
         MouseType = str_split(Sample, '_', simplify=TRUE)[,5],
         mouse_subj = paste0(Mouse,'_',MouseType,'_',Community)) %>%
  mutate(Diet = ifelse(Diet %in% 'noFiber','MD','Normal')) %>%
  mutate(Community = ifelse(Community %in% 'SC','StoolCom','DualCom')) %>%
  mutate(diet_com = paste0(Diet,'_',Community)) %>%
  mutate(diet_com = factor(diet_com, levels=c('Normal_StoolCom','MD_StoolCom','Normal_DualCom','MD_DualCom'))) 


alpha.div <- estimate_richness(ps.shannon)
sample_data(ps.shannon)$Shannon <- alpha.div$Shannon

df.shannon <- data.frame(ps.shannon@sam_data) 

stat_test <- df.shannon %>%
  group_by(Location) %>%
  wilcox_test(Shannon ~ diet_com) %>%
  adjust_pvalue(method="BH") %>%
  add_significance() %>%
  #filter(!p.adj.signif=='ns') %>%
  mutate(Community1 = str_split(group1,'_',simplify=TRUE)[,2],
         Community2 = str_split(group2,'_',simplify=TRUE)[,2]) %>%
  mutate(Diet1 = str_split(group1,'_',simplify=TRUE)[,1],
         Diet2 = str_split(group2,'_',simplify=TRUE)[,1]) %>%
  filter(Community1==Community2 | Diet1==Diet2) %>%
  add_xy_position(x="diet_com")
stat_test

# Get medians
df.shannon %>% group_by(Location, diet_com) %>%
  dplyr::summarise(shannon_median = median(Shannon))
# Plot
(g<-ggplot(df.shannon, aes(y=Shannon, x=diet_com, color=diet_com)) +
    geom_beeswarm()+
    geom_boxplot(alpha=0)+
    facet_wrap(~Location) +
    scale_color_manual(values=c('red','red','purple','purple')) +
    stat_pvalue_manual(stat_test, label = "p.adj.signif",
                       tip.length = 0.01,
                       bracket.shorten = 0.05) +
    ylim(0,5) +
    ylab('Shannon diversity for each mouse') +
    scale_x_discrete(labels=c('SD','MD','SD','MD'))+
    theme_minimal ()  )
ggsave(paste0(fig_dir, 'subpanels/Fig_2F_shannonDiv.pdf'),width=4,height=3)


#######################################
# Fig. 2G: Ordination analysis to show similarity analysis of StoolCom and DualCom colonization states
#######################################

mouse.on.ps <- ps %>% subset_samples(Project %in% 'mouse_exp2' & !SampleType %in% c('Cecum','Control'))
sample_data(mouse.on.ps) <- data.frame(mouse.on.ps@sam_data) %>%
  mutate(rep_merged= paste0(Location,'_',Diet,'_',Community))
mouse.on.ps <- mouse.on.ps %>% merge_samples("rep_merged") 

### Correlation analysis between two mouse colonization conditions
ps.corr <-mouse.on.ps %>% 
  prune_taxa(taxa_sums(.) > 0, .) %>%
  transform_sample_counts( function(x)  log10(x + 0.001))

asvs2correlate <- otu_table(ps.corr)
rownames(asvs2correlate) <- rownames(sample_data(ps.corr))

res <- rcorr(t(asvs2correlate), type="pearson") # compute correlations
corr.df.tmp <- flattenCorrMatrix(res$r, res$P) 
p.vals <- p.adjust(corr.df.tmp$p ,method="bonferroni") # adjust pvalues
r.vals <- corr.df.tmp$cor
corr.df <- cbind(corr.df.tmp %>% select(-p), p.vals)
corr.df

### Ordination
mouse.ps <- ps %>% subset_samples(Project %in% 'mouse_exp2' & Location %in% c('Stool','Small Intestine')) %>%
  merge_samples("mouse_location")

Sample<-rownames(mouse.ps@sam_data)

sample_data(mouse.ps) <- data.frame(mouse.ps@sam_data) %>% 
  cbind(Sample) %>%
  mutate(Mouse = str_split(Sample, '_', simplify=TRUE)[,1],
         Location = str_split(Sample, '_', simplify=TRUE)[,2],
         Diet = str_split(Sample, '_', simplify=TRUE)[,3],
         Community = str_split(Sample, '_', simplify=TRUE)[,4],
         MouseType = str_split(Sample, '_', simplify=TRUE)[,5],
         mouse_subj = paste0(Mouse,'_',MouseType,'_',Community),
         com_loc = paste0(Community,'_',Location)) %>%
  mutate(Diet = ifelse(Diet %in% 'noFiber','MD','Normal')) %>%
  mutate(Diet = factor(Diet, levels=c('Normal','MD'))) 

taxa2keep <- mouse.ps %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% filter(Abundance > 0.001)

mouse.ord.ps <- prune_taxa(taxa2keep$OTU, mouse.ps)  %>% transform_sample_counts(function(x) x/sum(x))

pcoa.bray <- ordinate(mouse.ord.ps,  method = "MDS", distance = "bray")

var_exp <- get_evals(pcoa.bray)$variance_exp

plot_ordination(mouse.ord.ps, pcoa.bray, shape="Diet", color="Community") +
    geom_point(size=1)+
    scale_color_manual(values=c('red','purple'))+
    scale_shape_manual(values=c(1,19)) +
    theme_minimal()+
    coord_fixed(sqrt(var_exp[2] / var_exp[1])) +
    facet_wrap(~Location, nrow=2)
ggsave(paste0(fig_dir,'subpanels/Fig_2G_bray_ordination.pdf'), width=7, height = 4)

#######################################
# Analysis comparing mouse intestinal microbiota to human intestinal microbiota
#######################################

# Figure out which taxa are enriched in the small intestine compared to the stool of these mice
mouse.ps <- ps %>% subset_samples(Project %in% 'mouse_exp2' & Location %in% c('Stool','Small Intestine')) %>%
  subset_samples(Diet %in% 'Normal' & Community %in% 'SC') %>% #'Stool&SI'
  merge_samples("mouse_location")
Sample<-rownames(mouse.ps@sam_data)
sample_data(mouse.ps) <- data.frame(mouse.ps@sam_data) %>% 
  cbind(Sample) %>%
  mutate(Mouse = str_split(Sample, '_', simplify=TRUE)[,1],
         Location = str_split(Sample, '_', simplify=TRUE)[,2],
         Diet = str_split(Sample, '_', simplify=TRUE)[,3],
         Community = str_split(Sample, '_', simplify=TRUE)[,4],
         MouseType = str_split(Sample, '_', simplify=TRUE)[,5],
         mouse_subj = paste0(Mouse,'_',MouseType,'_',Community),
         com_loc = paste0(Community,'_',Location)) %>%
  mutate(Diet = ifelse(Diet %in% 'noFiber','MD','Normal')) %>%
  mutate(Diet = factor(Diet, levels=c('Normal','MD'))) 
# Filter for taxa that are detected in 0.1% abundance in at least 2 mice
taxa2keep <- mouse.ps %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% filter(Abundance > 0.001) %>% group_by(OTU) %>%
  dplyr::summarise(n=n()) %>% filter(n>1)

ps_da <- prune_taxa(taxa2keep$OTU, mouse.ps)

## Estimate size factors
dds <- phyloseq_to_deseq2(ps_da, design = ~ 1)
dds <- estimateSizeFactors(dds, type = "poscount")
size_factors <- sizeFactors(dds)
summary(size_factors)

# Since we do not want to artificially inflate the counts in voom transformation
# we normalize the size factors so the the smallest is equal to 1.
size_factors <- size_factors/min(size_factors)
summary(size_factors)
pstips<-ps_da
seqtab <- t(as(pstips@otu_table, "matrix"))
smpdata <- data.frame(pstips@sam_data) %>%
  mutate(Location = factor(Location, levels = c("Stool", "Small Intestine")))

taxtable <- data.frame(pstips@tax_table)  %>%
  # dplyr::select(Kingdom:Species) %>%
  rownames_to_column("SeqName") %>%
  mutate(OrgName = SeqName)

dsgn <- ~ Location

limmaRes <- limma_fit(
  seqtab, smpdata, dsgn, 
  sizefac = size_factors[rownames(smpdata)], 
  block = smpdata$mouse_subj, 
  taxtable = taxtable,
  alpha = 1)  %>%
  left_join(data.frame(pstips@tax_table) %>%
              rownames_to_column("SeqName"))

alpha <- 0.05
limmaRes_sig <- limmaRes %>%
  filter(adj.P.Val < alpha, abs(logFC) > 0.75) %>%
  mutate(OrgName = paste0(ifelse(is.na(Genus), Family, Genus), " ", Species, " (", SeqName, ")"),
         col = factor(ifelse(logFC > 0.75, "Increased in SI", "Increased in stool"))) %>%
  select(OrgName, logFC, adj.P.Val,col) %>% arrange(col)
limmaRes_sig %>% filter(col %in% 'Increased in SI')

require("ggrepel")

limmaRes <- limmaRes %>%
  mutate(Genus = ifelse(is.na(Genus)==TRUE, Family,Genus)) %>%
  mutate(Species = ifelse(is.na(Species)==TRUE, 'sp', Species)) %>%
  left_join(genus.order, by='Family') %>%
  mutate(colors = ifelse(adj.P.Val > 0.05, NA,colors)) %>%
  mutate(colors = ifelse(abs(logFC) < 0.75, NA,colors)) %>%
  mutate(OrgName = ifelse(adj.P.Val > 0.05 | abs(logFC) < 0.75, '', Genus)) 


ggplot(limmaRes, aes(y=-log10(adj.P.Val), x=logFC, label=OrgName)) +
  geom_point(fill=limmaRes$colors, size=3,pch=21, color='black')+
  geom_hline(yintercept=1.3, linetype="dashed")+
  geom_vline(xintercept=0.75, linetype="dashed")+
  geom_vline(xintercept=-0.75, linetype="dashed")+
  #ylim(0,9) +
  #xlim(-4,4) +
  #geom_text_repel()+
  theme_minimal() +
  theme(axis.text.x= element_text(color='black'),axis.text.y= element_text(color='black'))

#######################################
# Analysis for differentially enriched taxa during a diet switch
#######################################

# Figure out which taxa are enriched in the small intestine compared to the stool of these mice
mouse.ps <- ps %>% subset_samples(Project %in% 'mouse_exp2' & Location %in% c('Small Intestine')) %>% #,'Small Intestine'
  subset_samples(Community %in% 'Stool&SI') %>% 
  #subset_samples(Community %in% 'SC') %>%
  merge_samples("mouse_location")
Sample<-rownames(mouse.ps@sam_data)
sample_data(mouse.ps) <- data.frame(mouse.ps@sam_data) %>% 
  cbind(Sample) %>%
  mutate(Mouse = str_split(Sample, '_', simplify=TRUE)[,1],
         Location = str_split(Sample, '_', simplify=TRUE)[,2],
         Diet = str_split(Sample, '_', simplify=TRUE)[,3],
         Community = str_split(Sample, '_', simplify=TRUE)[,4],
         MouseType = str_split(Sample, '_', simplify=TRUE)[,5],
         mouse_subj = paste0(Mouse,'_',MouseType,'_',Community),
         com_loc = paste0(Community,'_',Location)) %>%
  mutate(Diet = ifelse(Diet %in% 'noFiber','MD','Normal')) %>%
  mutate(Diet = factor(Diet, levels=c('Normal','MD'))) 
# Filter for taxa that are detected in 0.1% abundance in at least 2 mice
taxa2keep <- mouse.ps %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% filter(Abundance > 0.001) %>% group_by(OTU) %>%
  dplyr::summarise(n=n()) %>% filter(n>1)

ps_da <- prune_taxa(taxa2keep$OTU, mouse.ps)

## Estimate size factors
dds <- phyloseq_to_deseq2(ps_da, design = ~ 1)
dds <- estimateSizeFactors(dds, type = "poscount")
size_factors <- sizeFactors(dds)
summary(size_factors)

# Since we do not want to artificially inflate the counts in voom transformation
# we normalize the size factors so the the smallest is equal to 1.
size_factors <- size_factors/min(size_factors)
summary(size_factors)
pstips<-ps_da
seqtab <- t(as(pstips@otu_table, "matrix"))
smpdata <- data.frame(pstips@sam_data) 
#%>%
#  mutate(Location = factor(Location, levels = c("Stool", "Small Intestine")))

taxtable <- data.frame(pstips@tax_table)  %>%
  # dplyr::select(Kingdom:Species) %>%
  rownames_to_column("SeqName") %>%
  mutate(OrgName = SeqName)

dsgn <- ~ Diet

limmaRes <- limma_fit(
  seqtab, smpdata, dsgn, 
  sizefac = size_factors[rownames(smpdata)], 
  block = smpdata$mouse_subj, 
  taxtable = taxtable,
  alpha = 1)  %>%
  left_join(data.frame(pstips@tax_table) %>%
              rownames_to_column("SeqName"))

alpha <- 0.05
limmaRes_sig <- limmaRes %>%
  filter(adj.P.Val < alpha, abs(logFC) > 0.75) %>%
  mutate(OrgName = paste0(ifelse(is.na(Genus), Family, Genus), " ", Species, " (", SeqName, ")"),
         col = factor(ifelse(logFC > 0.75, "Increased in MD", "Increased in normal"))) %>%
  select(OrgName, logFC, adj.P.Val,col) %>% arrange(col)
limmaRes_sig
#write.csv(limmaRes_sig, paste0(clean_data_dir,'DualCom_Stool_differentiallyEnrichedAcrossDiet.csv'))
require("ggrepel")

limmaRes <- limmaRes %>%
  mutate(Genus = ifelse(is.na(Genus)==TRUE, Family,Genus)) %>%
  mutate(Species = ifelse(is.na(Species)==TRUE, 'sp', Species)) %>%
  left_join(genus.order, by='Family') %>%
  mutate(colors = ifelse(adj.P.Val > 0.05, NA,colors)) %>%
  mutate(colors = ifelse(abs(logFC) < 0.75, NA,colors)) %>%
  mutate(OrgName = ifelse(adj.P.Val > 0.05 | abs(logFC) < 0.75, '', Genus)) 


ggplot(limmaRes, aes(y=-log10(adj.P.Val), x=logFC, label=OrgName)) +
  geom_point(fill=limmaRes$colors, size=3,pch=21, color='black')+
  geom_hline(yintercept=1.3, linetype="dashed")+
  geom_vline(xintercept=0.75, linetype="dashed")+
  geom_vline(xintercept=-0.75, linetype="dashed")+
  #ylim(0,9) +
  #xlim(-4,4) +
  geom_text_repel()+
  theme_minimal() +
  theme(axis.text.x= element_text(color='black'),axis.text.y= element_text(color='black'))


#######################################
# Fig. 2H: Enterococcus ASVs aact differently in different conditions
#######################################

mouse.ps <- ps %>% subset_samples(Project %in% 'mouse_exp2' & Location %in% c('Stool','Small Intestine')) %>%
  merge_samples("mouse_location")

Sample<-rownames(mouse.ps@sam_data)

sample_data(mouse.ps) <- data.frame(mouse.ps@sam_data) %>% 
  cbind(Sample) %>%
  mutate(Mouse = str_split(Sample, '_', simplify=TRUE)[,1],
         Location = str_split(Sample, '_', simplify=TRUE)[,2],
         Diet = str_split(Sample, '_', simplify=TRUE)[,3],
         Community = str_split(Sample, '_', simplify=TRUE)[,4],
         MouseType = str_split(Sample, '_', simplify=TRUE)[,5],
         mouse_subj = paste0(Mouse,'_',MouseType,'_',Community),
         com_loc = paste0(Community,'_',Location)) %>%
  mutate(Diet = ifelse(Diet %in% 'noFiber','MD','Normal')) %>%
  mutate(Diet = factor(Diet, levels=c('Normal','MD'))) 

taxa2keep <- mouse.ps %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% filter(Abundance > 0.001)

mouse.df.entero <- prune_taxa(taxa2keep$OTU, mouse.ps)  %>% transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>%
  filter(Location %in% 'Small Intestine') %>%
  mutate(Community = ifelse(Community %in% 'SC','StoolCom', 'DualCom')) %>%
  mutate(Community = factor(Community, levels=c('StoolCom','DualCom'))) %>%
  mutate(Species = ifelse(is.na(Species)==TRUE, 'sp.', Species)) %>%
  mutate(Species_name = paste0(Genus,' ',Species,' ',OTU)) %>%
  mutate(Abundance = log10(Abundance + 0.001)) %>%
  filter(Genus %in% 'Enterococcus')

stat_test <- mouse.df.entero %>%
  group_by(Species_name, Community) %>%
  wilcox_test(Abundance ~ Diet) %>%
  adjust_pvalue(method="BH") %>%
  add_significance() %>%
  #filter(!p.adj.signif=='ns') %>%
  add_xy_position(x="Diet")
stat_test

ggplot(mouse.df.entero, aes(x=Diet, y=Abundance)) +
  geom_point() +
  geom_boxplot(alpha=0) +
  theme_minimal() +
  stat_pvalue_manual(stat_test, label = "p.adj.signif",
                     tip.length = 0.01,
                     bracket.shorten = 0.05) +
  facet_wrap(~Species_name+Community, nrow=3)

ggsave(paste0(fig_dir,'subpanels/Fig_2H_enteroAbundance.pdf'), width=3, height=4)
