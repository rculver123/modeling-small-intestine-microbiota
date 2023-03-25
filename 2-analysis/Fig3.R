#######################################
# The following script includes analysis for the repeatability of
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
# Fig. 3A: Schematic in illustrator of mice experiment
#######################################

#######################################
# Fig. S3A: Inoculum used in all mouse experiments was highly correlated within inoculum type
#######################################

ps2plot<- ps %>%
  subset_samples(Project %in% 'MouseExp2' | SampleType %in% c('Overnight')) 

taxa2keep<-psmelt(ps2plot %>% transform_sample_counts(function(x) x/sum(x))) %>%
  filter(Abundance > 0.001) 

ps2plot <- prune_taxa(taxa2keep$OTU, ps2plot) %>% transform_sample_counts(function(x) x/sum(x))

overnight.compare.df <- ps2plot %>% 
  psmelt() %>% 
  mutate(Family = factor(ifelse(is.na(Family), ifelse(!is.na(Phylum), paste0("(Unclassified ", Phylum, ")"), "(Unclassified phylum)"), Family))) %>%
  select(OTU, Sample, Abundance, Phylum, Family, Genus, Species) %>%
  #group_by(Sample, Family) %>%
  #dplyr::summarise(Abundance = sum(Abundance)) %>% 
  filter(Abundance > 0) %>%
  mutate(Abundance = log10(Abundance+0.001)) %>%
  select(OTU,Sample, Family, Abundance) %>%
  pivot_wider(values_from = Abundance, names_from = Sample, values_fill=-3) %>% data.frame() %>%
  left_join(genus.order, by="Family") 

colnames(overnight.compare.df)

# MouseExp2_p3_A8_Day0_Gavage__S200 was the inoculum to create StoolCom-colonized mice
# MouseExp2_p3_A7_Day0_Gavage__S199 was the SICom inoculum used in a 1:1 micture with the above StoolCom to create DualCom
# S318_mouse_exp1 is for the first gavage of StoolCom to create the StoolCom2-colonized mice
# S319_mouse_exp1 is for the first gavage of SICom to create the SICom-colonized mice
# S320_mouse_exp1 is the second gavage of StoolCom to create StoolCom->SICom-colonized mice
# S321_mouse_exp1 is the second gavage of SICom to create SICom->StoolCom2-colonized mice

# Compare StoolCom-colonized inoculum to StoolCom2-colonized inoculum
ggplot(overnight.compare.df,
           aes(x=MouseExp2_p3_A8_Day0_Gavage__S200, y=S318_mouse_exp1)) +
    geom_point(color=overnight.compare.df$colors) +
    geom_abline(intercept=0, slope=1) +
    theme_minimal() +
    scale_y_continuous(breaks=c(-1,-2,-3))+
    scale_x_continuous(breaks=c(-1,-2,-3))+
    stat_cor(method='pearson') +
    xlab('StoolCom used for StoolCom-colonized mice')+
    ylab('StoolCom used for StoolCom2-colonized mice')
ggsave(paste0(fig_dir,'subpanels/Fig_S3A_StoolComComparison.pdf'), width=2.5,height=2.5)

# Compare the SICom used to create DualCom-colonized mice and the SICom used to create SICom-mice
ggplot(overnight.compare.df,
       aes(x=MouseExp2_p3_A7_Day0_Gavage__S199, y=S319_mouse_exp1)) +
  geom_point(color=overnight.compare.df$colors) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  scale_y_continuous(breaks=c(-1,-2,-3))+
  scale_x_continuous(breaks=c(-1,-2,-3))+
  stat_cor(method='pearson') +
  xlab('SICom used for DualCom-colonized mice')+
  ylab('SICom used for second cross-colonization')
ggsave(paste0(fig_dir,'subpanels/Fig_S3B_SIComComparison.pdf'), width=2.5,height=2.5)

# Compare the two StoolCom's that were used to create SICom->StoolCom-colonized and StoolCom->SICom-colonized mice
ggplot(overnight.compare.df,
       aes(x=S320_mouse_exp1, y=S318_mouse_exp1)) +
  geom_point(color=overnight.compare.df$colors) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  scale_y_continuous(breaks=c(-1,-2,-3))+
  scale_x_continuous(breaks=c(-1,-2,-3))+
  stat_cor(method='pearson') +
  xlab('StoolCom used for StoolCom-colonized mice')+
  ylab('StoolCom used for second cross-colonization')
ggsave(paste0(fig_dir,'subpanels/Fig_S3C_StoolCom_crossComparison.pdf'), width=2.5,height=2.5)

# Compare the two SICom's that were used to create SICom->StoolCom-colonized and StoolCom->SICom-colonized mice
ggplot(overnight.compare.df,
       aes(x=S321_mouse_exp1, y=S319_mouse_exp1)) +
  geom_point(color=overnight.compare.df$colors) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  scale_y_continuous(breaks=c(-1,-2,-3))+
  scale_x_continuous(breaks=c(-1,-2,-3))+
  stat_cor(method='pearson') +
  xlab('SICom used for StoolCom-colonized mice')+
  ylab('SICom used for second cross-colonization')
ggsave(paste0(fig_dir,'subpanels/Fig_S3D_SICom_crossComparison.pdf'), width=2.5,height=2.5)



#######################################
# Supplementary Table 1: Stool and small intestine samples are highly correlated within a cage 
#######################################

mouse.ps <- ps %>%
  subset_samples(Project %in% c('mouse_exp1') & finalDaySamples %in% 1 & !SampleType %in% 'Cecum') %>%
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
  mutate(Community = ifelse(Community %in% 'SC', 'StoolCom',
                            ifelse(Community %in% 'IC', 'SICom',
                                   ifelse(Community %in% 'IC into SC-mouse', 'SICom->StoolCom-mouse', 'StoolCom ->SICom-mouse')))) %>%
  mutate(Diet = ifelse(Diet %in% 'noFiber','MD','SD')) %>%
  arrange(Location, Community, Diet)
write.csv(corr.df, paste0(fig_dir, 'Supplementary Table 2.csv'))

#######################################
# Fig. 3B: Barplots of SICom, StoolCom2, and cross-colonizations
#######################################

ps2plot<- ps %>% 
  subset_samples(Project %in% c('mouse_exp1') & finalDaySamples %in% 1 & !SampleType %in% 'Cecum') %>%
  merge_samples("mouse_location")

taxa2keep<-psmelt(ps2plot %>% transform_sample_counts(function(x) x/sum(x))) %>%
  filter(Abundance > 0.001) %>%
  group_by(OTU) %>% dplyr::summarise(num=n()) %>%
  filter(num>3) #must be in at least 2 mice

ps2plot <- prune_taxa(taxa2keep$OTU, ps2plot) %>% transform_sample_counts(function(x) x/sum(x))

df2plot <- psmelt(ps2plot) %>% 
  mutate(Mouse = str_split(Sample, '_', simplify=TRUE)[,1],
         Location = str_split(Sample, '_', simplify=TRUE)[,2],
         Diet = str_split(Sample, '_', simplify=TRUE)[,3],
         Community = str_split(Sample, '_', simplify=TRUE)[,4],
         MouseType = str_split(Sample, '_', simplify=TRUE)[,5],
         mouse_subj = paste0(Mouse,'_',MouseType,'_',Community),
         com_loc = paste0(Community,'_',Location)) %>%
  mutate(Community = factor(Community, levels=c('IC','SC into IC-mouse','SC','IC into SC-mouse'))) %>%
  mutate(Family = factor(ifelse(is.na(Family), ifelse(!is.na(Phylum), paste0("(Unclassified ", Phylum, ")"), "(Unclassified phylum)"), Family))) %>%
  left_join(genus.order, by="Family") %>%
  mutate(Family = factor(Family, levels=genus.order$Family)) %>%
  group_by(Sample, Location, Community, Family, colors) %>%
  dplyr::summarise(Abundance=sum(Abundance)) %>% ungroup() 

ggplot(df2plot , aes(x=Sample, y=Abundance, fill=Family)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=df2plot$colors) +
    theme_minimal() +
    theme(panel.margin.x=unit(0, "lines")) +
    theme(axis.text.y= element_text(color='black'), axis.text.x=element_blank(), axis.ticks.x=element_blank())+   
    ylab('Rel. ab.')   +           
    xlab('') +
    facet_wrap(~Location+Community, scale='free_x',nrow=2)                                                                      

ggsave(paste0(fig_dir, 'subpanels/Fig_3B_miceBarplots.pdf'), width = 10, height = 4)


# Priority effects visualization
ggplot(df2plot, aes(y=Family, x=Community, fill=log10(Abundance))) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~Location, scale='free_x')


#######################################
# Fig. 3C: Correlation between microbiota of mice in the small intestine and stool
# for StoolCom-colonized and StoolCom2-colonized mice
#######################################

# Get physeq object and merge all mice
sc.mice.ps <- ps %>% subset_samples(Project %in% 'mouse_exp2' & !SampleType %in% c('Cecum','Control') & Community %in% 'SC' |
                                       Project %in% c('mouse_exp1') & finalDaySamples %in% 1 & !SampleType %in% 'Cecum' & Community %in% 'SC')
sample_data(sc.mice.ps) <- data.frame(sc.mice.ps@sam_data) %>%
  mutate(rep_merged= paste0(Project,'_',Location,'_',Community))
sc.mice.ps <- sc.mice.ps %>% merge_samples("rep_merged") 

taxa2keep <- sc.mice.ps %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% filter(Abundance > 0.001)
sc.mice.ps<- prune_taxa(taxa2keep$OTU, sc.mice.ps)

rep.df <- sc.mice.ps %>% transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% 
  mutate(Family = factor(ifelse(is.na(Family), ifelse(!is.na(Phylum), paste0("(Unclassified ", Phylum, ")"), "(Unclassified phylum)"), Family))) %>%
  select(OTU, Sample, Abundance, Phylum, Family, Genus, Species) %>%
  #group_by(Sample, Family) %>%
  #dplyr::summarise(Abundance = sum(Abundance)) %>% 
  filter(Abundance > 0) %>%
  mutate(Abundance = log10(Abundance+0.001)) %>%
  select(OTU,Sample, Family, Abundance) %>%
  pivot_wider(values_from = Abundance, names_from = Sample, values_fill=-3) %>% data.frame() %>%
  left_join(genus.order, by="Family") 

colnames(rep.df)

t <- rep.df %>% filter(!mouse_exp1_Stool_SC == -3 && !mouse_exp2_Stool_SC==-3)
(a<-ggplot(t,
          aes(x=mouse_exp1_Stool_SC, y=mouse_exp2_Stool_SC)) +
  geom_point(color=t$colors) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  scale_y_continuous(breaks=c(-1,-2,-3))+
  scale_x_continuous(breaks=c(-1,-2,-3))+
  stat_cor(method='pearson') +
  xlab('Stool of \nStoolCom-colonized mice')+
  ylab('Stool of \nStoolCom2-colonized mice'))
ggsave(paste0(fig_dir,'subpanels/Fig_3C_stoolComparisons_a.pdf'), width=2.5,height=2.5)

t <- rep.df %>% filter(!mouse_exp1_Small.Intestine_SC == -3 && !mouse_exp2_Small.Intestine_SC==-3)
(b<-ggplot(t,
          aes(x=mouse_exp1_Small.Intestine_SC, y=mouse_exp2_Small.Intestine_SC)) +
  geom_point(color=t$colors) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  coord_fixed(1/1) +
  scale_y_continuous(breaks=c(-1,-2,-3))+
  scale_x_continuous(breaks=c(-1,-2,-3))+
  stat_cor(method='pearson') +
  xlab('Small intestine of \nStoolCom-colonized mice')+
  ylab('Small intestine of \nStoolCom2-colonized mice'))

ggsave(paste0(fig_dir,'subpanels/Fig_3C_stoolComparisons_b.pdf'), width=2.5,height=2.5)



#######################################
# Emergent ASVs analysis in StoolCom2
#######################################

# Get ASVs/families in the in vitro communities

overnight.cultures <- ps %>% subset_samples(Project %in% 'MouseExp2' | SampleType %in% 'Overnight') %>% 
  merge_samples("Community")
taxa2keep <- overnight.cultures %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% filter(Abundance > 0.001)
stool.com <- prune_taxa(taxa2keep$OTU, overnight.cultures) %>% transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% filter(Sample %in% 'SC', Abundance > 0.001)
unique(stool.com$Family)

#Mouse samples dataframe
sc.mice.ps <- ps %>% subset_samples(Project %in% 'mouse_exp2' & !SampleType %in% c('Cecum','Control') & Community %in% 'SC' |
                                      Project %in% c('mouse_exp1') & finalDaySamples %in% 1 & !SampleType %in% 'Cecum' & Community %in% 'SC')
sample_data(sc.mice.ps) <- data.frame(sc.mice.ps@sam_data) %>%
  mutate(rep_merged= paste0(Project,'_',Location,'_',Community))
sc.mice.ps <- sc.mice.ps %>% merge_samples("rep_merged") 

taxa2keep <- sc.mice.ps %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% filter(Abundance > 0.001)

mouse.df<- prune_taxa(taxa2keep$OTU, sc.mice.ps) %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>%
  filter(Abundance > 0.001)

# StoolCom2 mice emergent ASVs
stoolcom.stool <- mouse.df %>% filter(Sample=='mouse_exp2_Stool_SC') %>%
  filter(!OTU %in% stool.com$OTU) 
stoolcom.stool %>% group_by(Sample) %>% dplyr::summarise(Abund=sum(Abundance)) %>% dplyr::summarise(mean(Abund))

stoolcom.stool <- mouse.df %>% filter(Sample=='mouse_exp2_Small Intestine_SC') %>%
  filter(!OTU %in% stool.com$OTU) 
stoolcom.stool %>% group_by(Sample) %>% dplyr::summarise(Abund=sum(Abundance)) %>% dplyr::summarise(mean(Abund))

# What ASVs were detected in SI of StoolCom-colonized mice and not in StoolCom2?
stoolcom.only <- mouse.df %>% filter(Sample=='mouse_exp1_Small Intestine_SC') %>%
  filter(Abundance > 0.001) 
stoolcom2.only <- mouse.df %>% filter(Sample=='mouse_exp2_Small Intestine_SC') %>%
  filter(Abundance > 0.001) 

length(setdiff(stoolcom.only$OTU, stoolcom2.only$OTU))
stoolcom.asvs<-setdiff(stoolcom.only$OTU, stoolcom2.only$OTU)
stoolcom.only %>% filter(OTU %in% stoolcom.asvs) %>% group_by(Sample) %>% dplyr::summarise(Abund=sum(Abundance)) %>% dplyr::summarise(mean(Abund))
length(setdiff(stoolcom.only$Family, stoolcom2.only$Family))

length(setdiff(stoolcom2.only$OTU, stoolcom.only$OTU))
stoolcom2.asvs<-setdiff(stoolcom2.only$OTU, stoolcom.only$OTU)
stoolcom2.only %>% filter(OTU %in% stoolcom2.asvs) %>% group_by(Sample) %>% dplyr::summarise(Abund=sum(Abundance)) %>% dplyr::summarise(mean(Abund))
length(setdiff(stoolcom2.only$Family, stoolcom.only$Family))


#######################################
# Fig. 3D: Correlations between cross-colonizations
#######################################

ps.mouseexp1 <- ps %>% 
  subset_samples(Project %in% c('mouse_exp1') & finalDaySamples %in% 1 & !SampleType %in% 'Cecum')
sample_data(ps.mouseexp1) <- data.frame(ps.mouseexp1@sam_data) %>%
  mutate(rep_merged= paste0(Project,'_',Location,'_',Community))
ps.mouseexp1 <- ps.mouseexp1 %>% merge_samples("rep_merged") 

taxa2keep <- ps.mouseexp1 %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% filter(Abundance > 0.001) 
ps.mouseexp1.filt<- prune_taxa(taxa2keep$OTU, ps.mouseexp1)

df2plot.corr <- ps.mouseexp1.filt %>% transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% 
  mutate(Family = factor(ifelse(is.na(Family), ifelse(!is.na(Phylum), paste0("(Unclassified ", Phylum, ")"), "(Unclassified phylum)"), Family))) %>%
  select(OTU, Sample, Abundance, Phylum, Family, Genus, Species) %>%
  #group_by(Sample, Family) %>%
  #dplyr::summarise(Abundance = sum(Abundance)) %>% 
  filter(Abundance > 0) %>%
  mutate(Abundance = log10(Abundance+0.001)) %>%
  select(OTU,Sample, Family, Abundance) %>%
  pivot_wider(values_from = Abundance, names_from = Sample, values_fill=-3) %>% data.frame() %>%
  left_join(genus.order, by="Family") 

colnames(df2plot.corr)

# Small intestine of cross-colonization
ggplot(df2plot.corr,
          aes(x=mouse_exp1_Small.Intestine_IC.into.SC.mouse, y=mouse_exp1_Small.Intestine_SC.into.IC.mouse)) +
  geom_point(color=df2plot.corr$colors) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  scale_y_continuous(breaks=c(-1,-2,-3))+
  scale_x_continuous(breaks=c(-1,-2,-3))+
  stat_cor(method='pearson') +
  xlab('Small intestine of \nSICom into StoolCom-colonized mice')+
  ylab('Small intestine of \nStoolCom into SICom-colonized mice')
ggsave(paste0(fig_dir,'subpanels/Fig_3D_crossScatterPlot_a.pdf.pdf'), width=2.5,height=2.5)


# Compute shared abundance
stool.shared<-ps.mouseexp1.filt %>% transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% 
  filter(Sample %in% c('mouse_exp1_Stool_SC into IC-mouse','mouse_exp1_Stool_IC into SC-mouse')) %>%
  filter(Abundance > 0) %>%
  group_by(OTU) %>%
  dplyr::summarise(num = n()) %>% filter(num>1)
dim(stool.shared)
ps.mouseexp1.filt %>% transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% 
  filter(Sample %in% c('mouse_exp1_Stool_SC into IC-mouse','mouse_exp1_Stool_IC into SC-mouse')) %>%
  filter(Abundance > 0, !OTU %in% stool.shared$OTU) %>%
  group_by(Sample) %>% dplyr::summarise(sum(Abundance))


# Compute shared abundance
si.shared<-ps.mouseexp1.filt %>% transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% 
  filter(Sample %in% c('mouse_exp1_Small Intestine_SC into IC-mouse','mouse_exp1_Small Intestine_IC into SC-mouse')) %>%
  filter(Abundance > 0) %>%
  group_by(OTU) %>%
  dplyr::summarise(num = n()) %>% filter(num>1)
dim(si.shared)
ps.mouseexp1.filt %>% transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% 
  filter(Sample %in% c('mouse_exp1_Small Intestine_SC into IC-mouse','mouse_exp1_Small Intestine_IC into SC-mouse')) %>%
  filter(Abundance > 0, !OTU %in% stool.shared$OTU) %>%
  group_by(Sample) %>% dplyr::summarise(sum(Abundance))


# Stool of cross-colonization
ggplot(df2plot.corr,
       aes(x=mouse_exp1_Stool_IC.into.SC.mouse, y=mouse_exp1_Stool_SC.into.IC.mouse)) +
  geom_point(color=df2plot.corr$colors) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  scale_y_continuous(breaks=c(-1,-2,-3))+
  scale_x_continuous(breaks=c(-1,-2,-3))+
  stat_cor(method='pearson') +
  xlab('Stool of \nSICom into StoolCom-colonized mice')+
  ylab('Stool of \nStoolCom into ICCom-colonized mice')
ggsave(paste0(fig_dir,'subpanels/Fig_3D_crossScatterPlot_b.pdf.pdf'), width=2.5,height=2.5)


### SI only: Cross-colonizations comparison to DualCom

ps.mouseexp1 <- ps %>% 
  subset_samples(Project %in% c('mouse_exp1') & finalDaySamples %in% 1 & Location %in% 'Small Intestine' | 
                   Project %in% 'mouse_exp2' & Location %in% 'Small Intestine' & Diet %in% 'Normal' & Community %in% 'Stool&SI')
sample_data(ps.mouseexp1) <- data.frame(ps.mouseexp1@sam_data) %>%
  mutate(rep_merged= paste0(Project,'_',Location,'_',Community))
ps.mouseexp1 <- ps.mouseexp1 %>% merge_samples("rep_merged") 

taxa2keep <- ps.mouseexp1 %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% filter(Abundance > 0.001) 

ps.mouseexp1.filt<- prune_taxa(taxa2keep$OTU, ps.mouseexp1)

df2plot.corr <- ps.mouseexp1.filt %>% transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% 
  mutate(Family = factor(ifelse(is.na(Family), ifelse(!is.na(Phylum), paste0("(Unclassified ", Phylum, ")"), "(Unclassified phylum)"), Family))) %>%
  select(OTU, Sample, Abundance, Phylum, Family, Genus, Species) %>%
  #group_by(Sample, Family) %>%
  #dplyr::summarise(Abundance = sum(Abundance)) %>% 
  filter(Abundance > 0) %>%
  mutate(Abundance = log10(Abundance+0.001)) %>%
  select(OTU,Sample, Family, Abundance) %>%
  pivot_wider(values_from = Abundance, names_from = Sample, values_fill=-3) %>% data.frame() %>%
  left_join(genus.order, by="Family") 

colnames(df2plot.corr)

# Small intestine of cross-colonization
ggplot(df2plot.corr,
       aes(x=mouse_exp1_Small.Intestine_IC.into.SC.mouse, y=mouse_exp2_Small.Intestine_Stool.SI)) +
  geom_point(color=df2plot.corr$colors) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  scale_y_continuous(breaks=c(-1,-2,-3))+
  scale_x_continuous(breaks=c(-1,-2,-3))+
  stat_cor(method='pearson') +
  xlab('SI of SICom into \nStoolCom-colonized mice')+
  ylab('SI of DualCom- \ncolonized mice')
#ggsave(paste0(fig_dir,'subpanels/Fig_SX.pdf'), width=2.5,height=2.5)

ggplot(df2plot.corr,
       aes(x=mouse_exp1_Small.Intestine_SC.into.IC.mouse, y=mouse_exp2_Small.Intestine_Stool.SI)) +
  geom_point(color=df2plot.corr$colors) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  scale_y_continuous(breaks=c(-1,-2,-3))+
  scale_x_continuous(breaks=c(-1,-2,-3))+
  stat_cor(method='pearson') +
  xlab('SI of StoolCom into \nSICom-colonized mice')+
  ylab('SI of DualCom- \ncolonized mice')
#ggsave(paste0(fig_dir,'subpanels/Fig_SX.pdf'), width=2.5,height=2.5)


# Small intestine of cross-colonization
ggplot(df2plot.corr,
       aes(x=mouse_exp1_Small.Intestine_IC.into.SC.mouse, y=mouse_exp1_Small.Intestine_SC)) +
  geom_point(color=df2plot.corr$colors) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  scale_y_continuous(breaks=c(-1,-2,-3))+
  scale_x_continuous(breaks=c(-1,-2,-3))+
  stat_cor(method='pearson') +
  xlab('Small intestine of \nSICom into StoolCom-colonized mice')+
  ylab('Small intestine of \nStoolCom2-colonized mice')
# ggsave(paste0(fig_dir,'subpanels/Fig_3D_crossScatterPlot_a.pdf.pdf'), width=2.5,height=2.5)


#######################################
# Fig. 3E: Shannon diversity
#######################################

# Create data frame where each sample has corresponding Shannon div
ps.shannon<- ps %>% 
  subset_samples(Project %in% c('mouse_exp1') & finalDaySamples %in% 1 & !SampleType %in% 'Cecum') %>%
  merge_samples("mouse_location") %>%
  prune_taxa(taxa_sums(.) > 50, .) # because we are inserting information about raw reads, we do want to apply some sort of filter
# Since we filter for >5000 reads, we'll say taxa sums must be at least 50 (0.1%)
# Though without the filter, the data does not change_taxa()

alpha.div <- estimate_richness(ps.shannon)
sample_data(ps.shannon)$Shannon <- alpha.div$Shannon

df.shannon <- data.frame(ps.shannon@sam_data) %>%
  select(-Community, -Location) %>%
  rownames_to_column("Sample") %>%
  mutate(Community = str_split(Sample, '_',simplify=TRUE)[,4],
         Location = str_split(Sample,'_',simplify=TRUE)[,2]) %>%
  mutate(Community = factor(Community, levels=c('IC','SC into IC-mouse','SC','IC into SC-mouse'))) %>%
  select(Sample, Community, Location, Shannon) %>% data.frame()

# Number of mice for each category
df.shannon %>% group_by(Community, Location) %>% dplyr::summarise(n())

stat_test <- df.shannon %>%
  group_by(Location) %>%
  wilcox_test(Shannon ~ Community) %>%
  adjust_pvalue(method="BH") %>%
  add_significance() %>%
  filter(!p.adj.signif=='ns') %>%
  add_xy_position(x="Community")
stat_test

# Plot
(d<-ggplot(df.shannon , aes(y=Shannon, x=Community)) +
    geom_beeswarm(alpha=0.5) +
    geom_boxplot(alpha=0) +
    stat_pvalue_manual(stat_test, label = "p.adj.signif",
                       tip.length = 0.01,
                       bracket.shorten = 0.05) +
    ylab('Shannon diversity for each mouse') +
    theme_minimal () +
    ylim(0,4.5) +
    scale_x_discrete(labels=c('SICom','StoolCom into SICom','StoolCom','SICom into StoolCom'))+
    facet_wrap(~Location, nrow=1))

ggsave(paste0(fig_dir, 'subpanels/Fig_3E_shannonDiv.pdf'),width=3,height=2)

#######################################
# Fig. 3F: Ordination
#######################################

mouse.ps <- ps %>% 
  subset_samples(Project %in% c('mouse_exp1') & finalDaySamples %in% 1 & !SampleType %in% 'Cecum') %>%
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

(e<-plot_ordination(mouse.ord.ps, pcoa.bray, color="Community",shape="Community") +
    geom_point(size=1)+
    scale_color_manual(values=c('blue','plum','red','brown4'))+
    scale_shape_manual(values=c(1,1,1,1))+
    theme_minimal()+
    coord_fixed(sqrt(var_exp[2] / var_exp[1])) +
    facet_wrap(~Location, nrow=1))
ggsave(paste0(fig_dir,'subpanels/Fig_3F_bray_ordination.pdf'), width=7, height = 4)

#######################################
# Correlation with human intestine microbiota to the cross-colonization
#######################################

# Get physeq object and merge all mice
cross.mice.ps <- ps %>% subset_samples(Project %in% c('mouse_exp1') & finalDaySamples %in% 1 & !SampleType %in% 'Cecum' & !Community %in% c('SC','IC')|
                                      SampleType %in% 'original_2522_capsule')
sample_data(cross.mice.ps) <- data.frame(cross.mice.ps@sam_data) %>%
  mutate(rep_merged= paste0(Project,'_',Location,'_',Community))
cross.mice.ps <- cross.mice.ps %>% merge_samples("rep_merged") 

taxa2keep <- cross.mice.ps %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% filter(Abundance > 0.001)
cross.mice.ps<- prune_taxa(taxa2keep$OTU, cross.mice.ps)

rep.df <- cross.mice.ps %>% transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% 
  mutate(Family = factor(ifelse(is.na(Family), ifelse(!is.na(Phylum), paste0("(Unclassified ", Phylum, ")"), "(Unclassified phylum)"), Family))) %>%
  select(OTU, Sample, Abundance, Phylum, Family, Genus, Species) %>%
  group_by(Sample, Family) %>%
  dplyr::summarise(Abundance = sum(Abundance)) %>% 
  filter(Abundance > 0) %>%
  mutate(Abundance = log10(Abundance+0.001)) %>%
  select(Sample, Family, Abundance) %>%
  pivot_wider(values_from = Abundance, names_from = Sample, values_fill=-3) %>% data.frame() %>%
  left_join(genus.order, by="Family") 

colnames(rep.df)

t <- rep.df %>% filter(!mouse_exp1_Small.Intestine_SC.into.IC.mouse == -3 && !MouseExp3_NA_Human_Intestinal_Microbiota==-3)
(a<-ggplot(t,
           aes(x=mouse_exp1_Small.Intestine_SC.into.IC.mouse, y=MouseExp3_NA_Human_Intestinal_Microbiota)) +
    geom_point(color=t$colors) +
    geom_abline(intercept=0, slope=1) +
    theme_minimal() +
    scale_y_continuous(breaks=c(-1,-2,-3))+
    scale_x_continuous(breaks=c(-1,-2,-3))+
    stat_cor(method='pearson') +
    xlab('Small intestine of \nStoolCom into SICom-colonized mice')+
    ylab('Human small intestine 2522'))

t <- rep.df %>% filter(!mouse_exp1_Small.Intestine_IC.into.SC.mouse == -3 && !MouseExp3_NA_Human_Intestinal_Microbiota==-3)
(a<-ggplot(t,
           aes(x=mouse_exp1_Small.Intestine_IC.into.SC.mouse, y=MouseExp3_NA_Human_Intestinal_Microbiota)) +
    geom_point(color=t$colors) +
    geom_abline(intercept=0, slope=1) +
    theme_minimal() +
    scale_y_continuous(breaks=c(-1,-2,-3))+
    scale_x_continuous(breaks=c(-1,-2,-3))+
    stat_cor(method='pearson') +
    xlab('Small intestine of \nSICom into StoolCom-colonized mice')+
    ylab('Human small intestine 2522'))

ggsave(paste0(fig_dir, 'subpanels/Fig_2D_scatterPlots.pdf'), width = 5, height = 5)
