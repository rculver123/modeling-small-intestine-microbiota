#######################################
# The following script includes analysis for the
# passaging of small intestine and stool samples
# in media to generate in vitro communities
# 
# Rebecca Culver
#######################################

# configure directories, load libraries and base functions
rm(list=ls())
source(paste0(here::here(), "/0-config.R"))
ps<-readRDS(paste0(clean_data_dir,'merge_phyloseq.RDS'))
genus.order<-readRDS(paste0(fig_dir, 'official_colors.RDS'))

#######################################
# Fig. 1A Schematic of device sampling and collection
#######################################

#######################################
# Fig. 1B Barplots of initial inoculum and final day samples of in vitro communities
#######################################

ps2plot<- ps %>% 
  subset_samples(Project %in% 'c2m5' | SampleType %in% 'Human_Stool') %>% 
  subset_samples(Inoculum %in% c('1062','1064','1065','197') & Day %in% c('d0','d14') & !pH %in% c('5')| SampleType %in% 'Human_Stool')

taxa2keep<-psmelt(ps2plot %>% transform_sample_counts(function(x) x/sum(x))) %>%
  filter(Abundance > 0.001) %>% group_by(OTU) %>% 
  dplyr::summarise(num_OTU = n()) 

ps2plot <- prune_taxa(taxa2keep$OTU, ps2plot) %>% transform_sample_counts(function(x) x/sum(x))

df2plot <- psmelt(ps2plot) %>% 
  mutate(Inoculum = ifelse(is.na(Inoculum), 197, Inoculum)) %>%
  mutate(pH_SeqName = paste0(pH,'_',SeqName)) %>%
  mutate(Family = factor(ifelse(is.na(Family), ifelse(!is.na(Phylum), paste0("(Unclassified ", Phylum, ")"), "(Unclassified phylum)"), Family))) %>%
  left_join(genus.order, by="Family") %>%
  mutate(Family = factor(Family, levels=genus.order$Family)) %>%
  group_by(pH_SeqName, Inoculum, pH, Family, colors) %>%
  dplyr::summarise(Abundance=sum(Abundance)) %>% ungroup() %>%
  arrange(pH)

x_label <- df2plot$pH_SeqName %>% unique()

df2plot <- df2plot %>%
  mutate(pH_SeqName = factor(pH_SeqName, levels=x_label))

(b<-ggplot(df2plot , aes(x=pH_SeqName, y=Abundance, fill=Family)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=df2plot$colors) +
  theme_minimal() +
  theme(panel.margin.x=unit(0, "lines")) +
  theme(axis.text.y= element_text(color='black'), axis.text.x=element_blank(), axis.ticks.x=element_blank())+   
  ylab('Rel. ab.')   +           
  xlab('') +
  facet_wrap(~Inoculum, scale='free_x',nrow=1))                                                                       

ggsave(paste0(fig_dir, 'subpanels/Fig_1B_barplots.pdf'), width = 14, height = 4)

# Used to determine x axis labels
unique(df2plot %>% arrange(Inoculum) %>% select(pH_SeqName) %>% data.frame())

#######################################
# Fig. S1A: pH 5 generates communities dominated by a single species, often Lactobacillus
#######################################

pH5.ps <- ps %>% subset_samples (Project %in% 'c2m5' & Day %in% c('d14')) %>%
  subset_samples(pH %in% c(4,5)) 

taxa2keep<-psmelt(pH5.ps %>%
                    transform_sample_counts(function(x) x/sum(x))) %>%
  filter(Abundance > 0.001) %>% group_by(OTU) %>% 
  dplyr::summarise(num_OTU = n()) #%>%
 # filter(num_OTU > 2) # must be in at least 2 samples - REMOVING THIS PARAMETER DUE TO LACK OF REPEATABILITY IN THIS CONDITION

pH5.ps.filt <- prune_taxa(taxa2keep$OTU, pH5.ps)  %>% transform_sample_counts(function(x) x/sum(x))

# Barplots
df2plot <- psmelt(pH5.ps.filt) %>% 
  mutate(Family = factor(ifelse(is.na(Family), ifelse(!is.na(Phylum), paste0("(Unclassified ", Phylum, ")"), "(Unclassified phylum)"), Family))) %>%
  left_join(genus.order, by="Family") %>%
  mutate(Family = factor(Family, levels=genus.order$Family)) %>%
  group_by(SeqName, Inoculum, pH, Family, colors, Day) %>%
  dplyr::summarise(Abundance=sum(Abundance))

# Samples (n's)
dim(taxa2keep) # number of OTUs

(sa<-ggplot(df2plot , aes(x=SeqName, y=Abundance, fill=Family)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values= df2plot$colors) +
  theme_minimal() +
  theme(panel.margin.x=unit(0, "lines")) +
  theme(axis.text.y= element_text(color='black'), axis.text.x=element_blank(), axis.ticks.x=element_blank())+                                                                                                                                            
  ylab('Rel. ab.')   +                                                                                                                                                                                                                                                                                         
  facet_wrap(~Inoculum + pH, scale='free_x',nrow=1))     

ggsave(paste0(fig_dir, 'subpanels/Fig_S1A_pH5.pdf'), width = 5, height = 3)


#######################################
# Text: Final OD600 across communities after 48 hrs of growth
#######################################

## BECCA FILL THIS OUT HERE


#######################################
# Fig. S1B: Community composition equilibriated quickly
#######################################

# Stability: We have days 2, 4, 6, 8, and 14 (14 is the final day)
c2m5.ps <- ps %>% subset_samples (Project %in% 'c2m5' & Day %in% c('d2','d4','d6','d8','d14')) %>%
  subset_samples(!pH %in% c(0,4,5)) 

taxa2keep<-psmelt(c2m5.ps %>%
                    transform_sample_counts(function(x) x/sum(x))) %>%
  filter(Abundance > 0.001) %>% group_by(OTU) %>% 
  dplyr::summarise(num_OTU = n()) %>%
  filter(num_OTU > 1) # must be in at least 2 samples

c2m5.ps.filt <- prune_taxa(taxa2keep$OTU, c2m5.ps)  

c2m5.samdf <- sample_data(data.frame(c2m5.ps.filt@sam_data) %>% 
  mutate(Day = factor(Day, levels=c('d2','d4','d6','d8','d14'))) %>%
  mutate(Sample = paste0(Inoculum,"_", pH)))

# Samples (n's)
dim(taxa2keep) # number of OTUs
table(c2m5.samdf$pH) # sample number by pH
table(c2m5.samdf$Inoculum) # sample number by inoculum

# Reassign sample dataframe
c2m5.ps.filt@sam_data <- c2m5.samdf

# Barplots
df2plot <- psmelt(c2m5.ps.filt %>%
                    transform_sample_counts(function(x) x/sum(x))) %>% 
  mutate(Family = factor(ifelse(is.na(Family), ifelse(!is.na(Phylum), paste0("(Unclassified ", Phylum, ")"), "(Unclassified phylum)"), Family))) %>%
  left_join(genus.order, by="Family") %>%
  mutate(Family = factor(Family, levels=genus.order$Family)) %>%
  group_by(SeqName, Inoculum, pH, Family, colors, Day) %>%
  dplyr::summarise(Abundance=sum(Abundance)) %>%
  mutate(Inoculum = factor(Inoculum, levels=c(1062, 1064, 1065, 197))) %>%
  mutate(seq_ph=paste0(pH,'_',str_split(SeqName ,'_S',simplify=TRUE)[,2])) %>%
  mutate(pH = as.numeric(pH))


p<-ggplot(df2plot , aes(x=reorder(seq_ph,pH), y=Abundance, fill=Family)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=df2plot$colors) +
  theme_minimal() +
  theme(panel.margin.x=unit(0, "lines"), legend.position="bottom") +
  #theme(axis.text.y= element_text(color='black'), axis.text.x=element_blank(), axis.ticks.x=element_blank())+   
  ylab('Rel. ab.')   +   
  xlab('') +
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~Inoculum + Day, scale='free_x',ncol=5)     

ggsave(paste0(fig_dir, 'subpanels/Fig_S1B_stability.pdf'), p,width = 7, height = 7)

#######################################
# Fig. S1C and S1D: Community composition is repeatable across replicates
#######################################

# The following code calculates to pearson correlation of log10(relative abundance of OTUs)
# for all replicates

c2m5.ps.filt.noZeros <- c2m5.ps.filt %>%
  subset_samples(Day %in% 'd14') %>% #only looking at the last day
  prune_taxa(taxa_sums(.) > 0, .) %>% # remove any ASVs that are 0 now 
  transform_sample_counts(function(x) x/sum(x)) %>%
  transform_sample_counts(fun=filterfun1) %>% # turns all 0s into 0.001
  transform_sample_counts( function(x)  log10(x))

asvs2correlate <- otu_table(c2m5.ps.filt.noZeros)
rownames(asvs2correlate) <- sample_data(c2m5.ps.filt.noZeros)$Condition
  
res <- rcorr(t(asvs2correlate), type="pearson") # compute correlations
corr.df.tmp <- flattenCorrMatrix(res$r, res$P) 

p.vals <- p.adjust(corr.df.tmp$p ,method="bonferroni") # adjust pvalues
r.vals <- corr.df.tmp$cor

corr.df <- cbind(corr.df.tmp %>% select(-p), p.vals) %>% #rearraying dataframe to be easier to plot
  dplyr::rename(Sample_1=row, Sample_2=column, R=cor) %>%
  mutate(Inoculum_1 = str_split(Sample_1, '_', simplify=TRUE)[,2]) %>%
  mutate(pH_1 = str_split(Sample_1, '_', simplify=TRUE)[,1]) %>%
  mutate(Inoculum_2 = str_split(Sample_2, '_', simplify=TRUE)[,2]) %>%
  mutate(pH_2 = str_split(Sample_2, '_', simplify=TRUE)[,1]) %>%
  filter(Inoculum_1 == Inoculum_2, pH_1==pH_2, !R==1.000000) %>% # Remove comparisons with self
  mutate(pH_1 = factor(pH_1, levels=c(5,6,7,8,9,10))) %>%
  unique() %>% #each value is in there twice, so only take unique correlations
  mutate(sig = ifelse(p.vals < 0.001, 'black', 'red')) # determine sig correlations

table(corr.df$sig)

(sc<-ggplot(corr.df, aes(x= pH_1, y=R)) +
  geom_boxplot(alpha=0)+
  geom_beeswarm(color=corr.df$sig) +
  ylim(0,1) +
  theme_minimal() +
  ylab('Pearson value (R)') +
  xlab('pH of in vitro community') +
  facet_wrap(~Inoculum_1))
ggsave(paste0(fig_dir, 'subpanels/Fig_S1C_repeatability.pdf'), width = 4, height = 4)

### Example of a correlation for one device and one stool in vitro community (bonus info - not in paper)
                                                                   
sub.corr.df <- psmelt(c2m5.ps.filt.noZeros) %>%
  filter(Inoculum %in% c(1062,197), pH %in% 7) %>%
  select(OTU, Sample, Abundance, Family) %>%
  pivot_wider(names_from=Sample, values_from=Abundance) %>%
  left_join(genus.order, by="Family")
colnames(sub.corr.df)
  
(sd<-ggplot(sub.corr.df, aes(x=c2m5_ps2_S301, y=c2m5_ps2_S302)) +
  geom_point(color=sub.corr.df$colors) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  coord_fixed(1/1) +
  stat_cor(method='pearson') +
  xlab('Device #1062 pH 7 Rep 1')+
  ylab('Device #1062 pH 7 Rep 2'))

(sd2<-ggplot(sub.corr.df, aes(x=c2m5_ps2_S119, y=c2m5_ps2_S120)) +
  geom_point(color=sub.corr.df$colors) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  coord_fixed(1/1) +
  stat_cor(method='pearson') +
  xlab('Stool #197 pH 7 Rep 1')+
  ylab('Stool #197 pH 7 Rep 2'))

#######################################
# Fig. 1C: Shannon diversity
#######################################

# Create data frame where each sample has corresponding Shannon div
ps.shannon1<- ps %>% 
  subset_samples(Project %in% 'c2m5' & Day %in% c('d0','d14') & !SampleType %in% 'regrow' & !pH %in% 5| 
                   SampleType %in% c('Human_Stool')) %>%
  prune_taxa(taxa_sums(.) > 50, .) # because we are inserting information about raw reads, we do want to apply some sort of filter
  # Since we filter for >5000 reads, we'll say taxa sums must be at least 50 (0.1%)
  # Though without the filter, the data does not change

alpha.div <- estimate_richness(ps.shannon1)
sample_data(ps.shannon1)$Shannon <- alpha.div$Shannon

df.shannon1 <- data.frame(ps.shannon1@sam_data) %>% 
  
  mutate(pH = factor(pH, levels = c('0','6','7','8','9','10'))) %>%
  group_by(Inoculum, pH) %>%
  dplyr::summarise(avg_shannon = mean(Shannon),med_shannon = median(Shannon)) %>%
  mutate(Inoculum = as.factor(Inoculum))

# Plot
(c<-ggplot(df.shannon1, aes(y=avg_shannon, x=pH, group=Inoculum)) +
  geom_line(aes(linetype=Inoculum)) +
  geom_point(size=1) +
  ylab('Mean Shannnon diversity of replicates') +
  theme_minimal () +
  ylim(0,4))
ggsave(paste0(fig_dir, 'subpanels/Fig_1C_shannonDiv.pdf'),width=3,height=3)

#######################################
# Bonus explanantion: How correlated are the invitro communities with the initial inoculum?
# Looking at the log10(rel ab) of shared ASVs/families between the inoculum and invitro communities
# After adjusting p-values using bonferroni, there are no significant correlations
#######################################

# Create data frame to compare each set of 3 replicates with the initial inoculum
ps.compare <- ps %>% 
  subset_samples(Project %in% 'c2m5' & Day %in% c('d0','d14') & !SampleType %in% 'regrow' & !pH %in% 5| 
                   SampleType %in% c('Human_Stool')) %>%
  merge_samples("Condition") # add all read counts of replicates together

taxa2keep <- psmelt(ps.compare %>% transform_sample_counts(function(x) x/sum(x))) %>%
  filter(Abundance > 0.001) 

## ALL CORRELATIONS LOOK AT THE AVERAGE ACROSS REPLICATES COMPARED TO THE INOCULUM
## AND ARE COMPUTED BASED OFF LOG10(TAXA)

# Look first at correlations at the species-level 
compare.corr.df <- prune_taxa(taxa2keep$OTU, ps.compare) %>%
  transform_sample_counts(function(x) log10(x/sum(x))) %>%
  psmelt() %>% select(-Condition) %>% mutate(Condition=Sample) %>% #column renaming
  mutate(Abundance = ifelse(Abundance < -3, NA, Abundance)) %>%
  mutate(Abundance = ifelse(Abundance == "-Inf", NA, Abundance)) %>%
  select(Condition, OTU, Abundance) %>%
  pivot_wider(values_from = Abundance, names_from = OTU, values_fill=NA) %>% data.frame() #fill with -3 as that is the lowest ASV
rownames(compare.corr.df) <- compare.corr.df$Condition
compare.corr.df <- compare.corr.df %>% select(-Condition) %>% t()

# Determine correlation
res <- rcorr(compare.corr.df, type="pearson")

# Create dataframe of results
corr.df2plot <- flattenCorrMatrix(res$r, res$P) 
pvalues <- p.adjust(corr.df2plot$p ,method="bonferroni")
corr.df2plot.filt <- cbind(corr.df2plot %>% dplyr::rename(p_unadj=p), pvalues ) %>%
  dplyr::rename(p=pvalues) %>%
  mutate(p = ifelse(p < 0.05, p, 'not sig')) %>%
  mutate(p = ifelse(p == 0, '< 10e-15', p)) %>%
  arrange(row, column) %>%
  select(row, column, cor, p, p_unadj) %>% filter(!p=="not sig") %>%
  # Only get those that correlated with the initial inoculum 
  mutate(pH_1 = str_split(row, "_",simplify=TRUE)[,1] )%>%
  mutate(pH_2 = str_split(column, "_",simplify=TRUE)[,1] )%>%
  mutate(Inoculum1 = str_split(row, "_",simplify=TRUE)[,2] ) %>%
  mutate(Inoculum2 = str_split(column, "_",simplify=TRUE)[,2] ) %>%
  filter(pH_1==0 | pH_2 ==0, Inoculum1 == Inoculum2) %>% arrange(cor, decreasing=TRUE) %>%
  mutate(Inoculum=Inoculum1, pH = pmax(pH_1, pH_2))
corr.df2plot.filt

### FAMILY-LEVEL

# Because certain ASVs are annotated as NA at the family level, we first want to rename those to their Order
ps.renamed <- prune_taxa(taxa2keep$OTU, ps.compare)
taxa.renamed <- data.frame(tax_table(ps.renamed)) %>%
  mutate(Family = factor(ifelse(is.na(Family), ifelse(!is.na(Phylum), paste0("(Unclassified ", Phylum, ")"), "(Unclassified phylum)"), Family)))
tax_table(ps.renamed) <- as.matrix(taxa.renamed)

compare.corr.df <- ps.renamed %>%
  tax_glom("Family") %>%
  transform_sample_counts(function(x) log10(x/sum(x))) %>%
  psmelt() %>% select(-Condition) %>% mutate(Condition=Sample) %>% #column renaming
  mutate(Abundance = ifelse(Abundance < -3, NA, Abundance)) %>%
  mutate(Abundance = ifelse(Abundance == "-Inf", NA, Abundance)) %>%
  select(Condition, Family, Abundance) %>%
  pivot_wider(values_from = Abundance, names_from = Family, values_fill=NA) %>% data.frame() #fill with -3 as that is the lowest ASV
rownames(compare.corr.df) <- compare.corr.df$Condition
compare.corr.df <- compare.corr.df %>% select(-Condition) %>% t()

# Determine correlation
res <- rcorr(compare.corr.df, type="pearson")

# Create dataframe of results after correcting pvalues
corr.df2plot <- flattenCorrMatrix(res$r, res$P) 
pvalues <- p.adjust(corr.df2plot$p ,method="bonferroni")
corr.df2plot.filt <- cbind(corr.df2plot %>% dplyr::rename(p_unadj=p), pvalues ) %>%
  dplyr::rename(p=pvalues) %>%
  mutate(p_label = ifelse(p<0.001, '****',
                          ifelse(p < 0.001, '***',
                            ifelse(p<0.01, '**', 
                              ifelse(p<0.1, '*', ''))))) %>%
  arrange(row, column) %>%
  select(row, column, cor, p,p_unadj, p_label) %>% 
  filter(p<0.05) %>%
  # Only get those that correlated with the initial inoculum 
  mutate(pH_1 = str_split(row, "_",simplify=TRUE)[,1] )%>%
  mutate(pH_2 = str_split(column, "_",simplify=TRUE)[,1] )%>%
  mutate(Inoculum1 = str_split(row, "_",simplify=TRUE)[,2] ) %>%
  mutate(Inoculum2 = str_split(column, "_",simplify=TRUE)[,2] ) %>%
  filter(pH_1==0 | pH_2 ==0, Inoculum1 == Inoculum2) %>%
  mutate(pH = pmax(pH_1, pH_2)) %>%
  mutate(pH = factor(pH, levels=c(6,7,8,9,10))) %>%
  mutate(Inoculum = Inoculum1) %>% arrange(cor, decreasing=TRUE)
corr.df2plot.filt  

#######################################
# Text explanantion: How much of the family-level abundance is shared between in vitro communities and the inoculum?
# Calculate shared families between the inoculum and in vitro community, then calculate the summed relative abundance
# of those shared families in the inoculum
#######################################

# Replicates are not merged in this analysis by summing the total read abundance and then calculating the relative abundance
# ASVs that are under 0.1% relative abundance are not considered

# Create data frame to compare each set of 3 replicates with the initial inoculum
ps.shared <- ps %>% 
  subset_samples(Project %in% 'c2m5' & Day %in% c('d0','d14') & !SampleType %in% 'regrow' & !pH %in% 5| 
                   SampleType %in% c('Human_Stool')) %>%
  merge_samples("Condition") # add all read counts of replicates together

taxa2keep <- psmelt(ps.shared %>% transform_sample_counts(function(x) x/sum(x))) %>%
  filter(Abundance > 0.001) 

# Because certain ASVs are annotated as NA at the family level, we first want to rename those to their Order
ps.renamed <- prune_taxa(taxa2keep$OTU, ps.shared)
taxa.renamed <- data.frame(tax_table(ps.renamed)) %>%
  mutate(Family = factor(ifelse(is.na(Family), ifelse(!is.na(Phylum), paste0("(Unclassified ", Phylum, ")"), "(Unclassified phylum)"), Family)))
tax_table(ps.renamed) <- as.matrix(taxa.renamed)

shared.df <- prune_taxa(taxa2keep$OTU, ps.renamed) %>% transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% filter(Abundance > 0.001) %>%
  mutate(pH = as.numeric(str_split(Sample,'_', simplify=TRUE)[,1]) )#pH gets re-factored during merge, so we need to adjust back to regular pH values

# Families in the initial inoculum
stool.families <- shared.df %>% filter(Sample %in% '0_197')
length(stool.families$Family %>% unique())
d1062.families <- shared.df %>% filter(Sample %in% '0_1062')
length(d1062.families$Family %>% unique())
d1064.families <- shared.df %>% filter(Sample %in% '0_1064')
length(d1064.families$Family %>% unique())
d1065.families <- shared.df %>% filter(Sample %in% '0_1065')
length(d1065.families$Family %>% unique())

print('Stool-derived in vitro communities')
# Number of families total
shared.df %>% filter(Inoculum %in% 197, !pH %in% 0) %>%
  select(pH, Family) %>% unique() %>% group_by(pH) %>% dplyr::summarise(n())
# Number of families that were also in inoculum
shared.df %>% filter(Inoculum %in% 197, !pH %in% 0) %>%
  filter(Family %in% stool.families$Family) %>%
  select(pH, Family) %>% unique() %>% group_by(pH) %>% dplyr::summarise(n())
# Representative abundance at the family-level
shared.df %>% filter(Inoculum %in% 197, !pH %in% 0) %>%
  filter(Family %in% stool.families$Family) %>%
  group_by(Inoculum, pH) %>%
  dplyr::summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)

print('SI-derived in vitro communities')
# Number of families
shared.df %>% filter(!Inoculum %in% 197, !pH %in% 0) %>%
  select(pH, Family) %>% unique() %>% group_by(pH) %>% dplyr::summarise(n())
# Number of families that were also in inoculum
shared.df %>% filter(!Inoculum %in% 197, !pH %in% 0) %>%
  filter(Family %in% c(d1062.families$Family,d1064.families$Family,d1065.families$Family)) %>%
  select(pH, Inoculum, Family) %>% unique() %>% group_by(Inoculum,pH) %>% dplyr::summarise(n=n()) %>% arrange(n)

print('1062')
shared.df %>% filter(Inoculum %in% 1062, !pH %in% 0) %>%
  filter(Family %in% d1062.families$Family) %>%
  group_by(Inoculum, pH) %>%
  dplyr::summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
print('1064')
shared.df %>% filter(Inoculum %in% 1064, !pH %in% 0) %>%
  filter(Family %in% d1064.families$Family) %>%
  group_by(Inoculum, pH) %>%
  dplyr::summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
print('1065')
shared.df %>% filter(Inoculum %in% 1065, !pH %in% 0) %>%
  filter(Family %in% d1065.families$Family) %>%
  group_by(Inoculum, pH) %>%
  dplyr::summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)


#######################################
# Fig. 1D: Ordination of the inoculum and in vitro communities
#######################################

# In this plot, we will look at each individual replicates and inoculua at the log2(read count)
# Log2 was chosen to normalize the data and simultaneously account for heteroscedasticity
# Log10(relative abundance) cannot be used in Bray-curtis analysis due to negative numbers
# Only taxa above 0.1% abundance in at least 1 sample will be considered
# A Bray-curtis ordination will be used to stratify based on unique ASVs and their abundances

# Because we'll want to plot data from experiments discussed later in the paper on the pcoa axes,
# we're going to create a phyloseq object of all data 

ps.ordination <- ps %>% 
  subset_samples(Project %in% 'c2m5' & Day %in% c('d0','d14') & !SampleType %in% 'regrow' & !pH %in% 5 | 
                   SampleType %in% c('Human_Stool','DAY0_1','original_2522_capsule'))


taxa2keep <- psmelt(ps.ordination %>% transform_sample_counts(function(x) x/sum(x))) %>%
  filter(Abundance > 0.001) 

ps.ordination.filt <- prune_taxa(taxa2keep$OTU, ps.ordination)  %>% transform_sample_counts(function(x) x/sum(x))

ps.ordination.filt@sam_data <- sample_data(data.frame(ps.ordination.filt@sam_data) %>%
                                             mutate(Inoculum = as.factor(Inoculum)) %>%
                                             mutate(Buffer = ifelse(is.na(Buffer)==TRUE,0,Buffer)) %>%
                                             mutate(Buffer = factor(Buffer, levels=c(100,50,0))))

pcoa.bray <- ordinate(ps.ordination.filt,  method = "MDS", distance = "bray")

var_exp <- get_evals(pcoa.bray)$variance_exp
scores <- get_scores(pcoa.bray, sample_data(ps.ordination.filt))

pH.cols <-colorRampPalette(c('red','blue'))(8)
pH.cols<-c('#000000',pH.cols)

plot_ordination(ps.ordination.filt, pcoa.bray, shape="Inoculum", color="pH") +
  geom_point(size=1)+
  scale_color_manual(values=pH.cols)+
  scale_shape_manual(values=c(15,4,16,1,7)) +
  theme_minimal()+
  coord_fixed(sqrt(var_exp[2] / var_exp[1])) +
  facet_wrap(~Buffer, nrow=2)
ggsave(paste0(fig_dir,'subpanels/Fig_1D_bray_ordination.pdf'), width=7, height = 4)

#######################################
# Fig. 1E: Beta diversity analysis
#######################################

ps.ordination <- ps %>% 
  subset_samples(Project %in% 'c2m5' & Day %in% c('d0','d14') & !SampleType %in% 'regrow' & !pH %in% 5 | 
                   SampleType %in% c('Human_Stool','DAY0_1','original_2522_capsule'))

taxa2keep <- psmelt(ps.ordination %>% transform_sample_counts(function(x) x/sum(x))) %>%
  filter(Abundance > 0.001) 

ps.ordination.filt <- prune_taxa(taxa2keep$OTU, ps.ordination)  %>% transform_sample_counts(function(x) x/sum(x))

ps.ordination.filt@sam_data <- sample_data(data.frame(ps.ordination.filt@sam_data) %>%
                                             mutate(Inoculum = as.factor(Inoculum)) %>%
                                             mutate(Buffer = ifelse(is.na(Buffer)==TRUE,0,Buffer)) %>%
                                             mutate(Buffer = factor(Buffer, levels=c(0,50,100))))
df.min <- data.frame(ps.ordination.filt@sam_data) %>%
  select(SeqName, Inoculum, pH, Buffer)

dist.bray <- phyloseq::distance(ps.ordination.filt, method="bray", type="samples")

# Assess beta diversity achieved within a single inoculum using a variety of conditions
# vs beta diversity achieved by using different inocula
dist.df <- as.matrix(dist.bray) %>%
  data.frame() %>% 
  rownames_to_column("Sample1") %>%
  pivot_longer(!Sample1, names_to="Sample2", values_to="dist") %>%
  filter(Sample1 != Sample2, Sample1 < Sample2) %>%
  left_join(df.min, by=c("Sample1" = "SeqName")) %>%
  left_join(df.min, by=c("Sample2" = "SeqName"), suffix = c("_1","_2")) %>%
  # ONLY LOOK AT DEVICE SAMPLES
  filter(!Inoculum_1=='197', !Inoculum_2=='197') %>% data.frame()

dist.within <- dist.df %>%
  filter(Inoculum_1==Inoculum_2) %>% # since we're comparing within
  filter(pH_1 %in% c(6,6.5,7,7.5,8,8.5,9), pH_2 %in% c(6,6.5,7,7.5,8,8.5,9), Buffer_1 == 0, Buffer_2 == 0) %>%
  # Get average across replicates by goruping 
  group_by(Inoculum_1, Inoculum_2, pH_1, Buffer_1, pH_2, Buffer_2) %>%
  dplyr::summarise(dist=mean(dist)) %>% data.frame() %>%
  # Ignore replicates comparison to each other
  filter(!pH_1==pH_2) %>%
  mutate(category='Within inoculum')

dist.across <- dist.df %>%
  filter(!Inoculum_1==Inoculum_2) %>% # since we're comparing across inocula now
  filter(pH_1 %in% c(6,6.5,7,7.5,8,8.5,9), pH_2 %in% c(6,6.5,7,7.5,8,8.5,9), Buffer_1 == 0, Buffer_2 == 0) %>%
  # Get average across replicates by goruping 
  group_by(Inoculum_1, Inoculum_2, pH_1, Buffer_1, pH_2, Buffer_2) %>%
  dplyr::summarise(dist=mean(dist)) %>% data.frame() %>%
  # Ignore replicates comparison to each other
  filter(!pH_1==pH_2) %>%
  mutate(category='Across inocula')

dist.long <- rbind(dist.within, dist.across)

stat_test <- dist.long %>%
  wilcox_test(dist ~ category) %>%
  adjust_pvalue(method="bonferroni") %>%
  add_significance() %>%
  add_xy_position(x="category")
stat_test %>% data.frame()

# Plot the distance faceted by buffer
ggplot(dist.long, aes(y=dist, x=category)) +
  geom_boxplot(outlier.shape=NA) +
  geom_beeswarm(alpha=0.5) +
  stat_pvalue_manual(stat_test, label = "p.adj.signif",
                     tip.length = 0.01,
                     bracket.shorten = 0.05) +
  ylim(0,1) +
  theme_minimal()
ggsave(paste0(fig_dir, 'subpanels/Fig_1E_distBoxplot.pdf'),width=3,height=3)


# BONUS INFO NOT IN PAPER
# Measure the beta diversity of communities from inoculum

dist.df.buffer <- dist.df %>%
  filter(Inoculum_1==2522, Inoculum_2==2522) %>%
  # Get average across replicates by goruping 
  group_by(pH_1, Buffer_1, pH_2, Buffer_2) %>%
  dplyr::summarise(dist=mean(dist)) %>% ungroup() %>%  data.frame() %>%
  # Ignore replicates comparison to each other
  filter(!pH_1==pH_2) %>%
  mutate(pH_1 = as.numeric(as.character(pH_1)), pH_2=as.numeric(as.character(pH_2)) ) %>%
  # Get comparisons to inoculum only
  mutate(pH_min = pmin(pH_1, pH_2)) %>%
  mutate(pH = pmax(pH_1, pH_2)) %>%
  filter(pH_min == 0)
dist.df.buffer

stat_test <- dist.df.buffer %>%
  wilcox_test(dist ~ Buffer_1) %>%
  adjust_pvalue(method="bonferroni") %>%
  add_significance() %>%
  add_xy_position(x="Buffer_1")
stat_test

# Plot the distance faceted by buffer
ggplot(dist.df.buffer, aes(y=dist, x=Buffer_1)) +
    geom_boxplot(outlier.shape=NA) +
    geom_beeswarm(alpha=0.5) +
    stat_pvalue_manual(stat_test, label = "p.adj.signif",
                       tip.length = 0.01,
                       bracket.shorten = 0.05) +
    ylim(0,1) +
    ylab('Bray-curtis distance between different pH in vitro communities')+
    theme_minimal()


#######################################
# Fig. 1F: Barplots of buffered communities
# Text analysis: includes correlation of inoculum with buffered communities
#######################################

ps2plot<- ps %>% 
  subset_samples( !pH %in% c(5,5.5) & SampleType %in% c('DAY0_1','original_2522_capsule')) 

taxa2keep<-psmelt(ps2plot %>% transform_sample_counts(function(x) x/sum(x))) %>%
  filter(Abundance > 0.001) 

ps2plot <- prune_taxa(taxa2keep$OTU, ps2plot) %>% transform_sample_counts(function(x) x/sum(x))

df2plot <- psmelt(ps2plot) %>% 
  mutate(Buffer = ifelse(is.na(Buffer)==TRUE, -1, Buffer)) %>% # relabel inoculum as "-1" buffer
  mutate(Family = factor(ifelse(is.na(Family), ifelse(!is.na(Phylum), paste0("(Unclassified ", Phylum, ")"), "(Unclassified phylum)"), Family))) %>%
  left_join(genus.order, by="Family") %>%
  mutate(Family = factor(Family, levels=genus.order$Family)) %>%
  mutate(pH_SeqName = paste0(pH,'_',SeqName)) %>%
  group_by(pH_SeqName, Inoculum, pH, Buffer, Family, colors) %>%
  dplyr::summarise(Abundance=sum(Abundance)) %>% ungroup() %>%
  arrange(Buffer, pH)

# Plotting order
df2plot %>% filter(Abundance > 0) %>% select(pH, Buffer, Family) %>% unique() %>% data.frame()

(f<-ggplot(df2plot , aes(x=pH_SeqName, y=Abundance, fill=Family)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=df2plot$colors) +
    theme_minimal() +
    theme(panel.margin.x=unit(0, "lines")) +
    theme(axis.text.y= element_text(color='black'), axis.text.x=element_blank(), axis.ticks.x=element_blank())+   
    ylab('Rel. ab.')   +           
    xlab('') +
    facet_wrap(~Buffer, scale='free_x',nrow=1))                                                                       

ggsave(paste0(fig_dir, 'subpanels/Fig_1F_bufferBarplots.pdf'), width = 10, height = 4)

# Used to determine x axis labels
unique(df2plot %>% arrange(Buffer) %>% select(pH_SeqName) %>% data.frame())
#######################################
# Fig. 1G: Shannon diversity
#######################################

# Create data frame where each sample has corresponding Shannon div
ps.shannon2<- ps %>% 
  subset_samples( !pH %in% c(5,5.5) & SampleType %in% c('DAY0_1','original_2522_capsule')) %>%
  prune_taxa(taxa_sums(.) > 50, .) # because we are inserting information about raw reads, we do want to apply some sort of filter
# Since we filter for >5000 reads, we'll say taxa sums must be at least 50 (0.1%)
# Though without the filter, the data does not change_taxa()

alpha.div <- estimate_richness(ps.shannon2)
sample_data(ps.shannon2)$Shannon <- alpha.div$Shannon

df.shannon2 <- data.frame(ps.shannon2@sam_data) %>% 
  mutate(pH = factor(pH, levels = c('0','6','6.5','7','7.5','8','8.5','9'))) %>%
  group_by(Inoculum, pH, Buffer) %>%
  dplyr::summarise(avg_shannon = mean(Shannon),med_shannon = median(Shannon)) %>%
  mutate(Buffer = as.factor(Buffer))

# Plot
(g<-ggplot(df.shannon2, aes(y=avg_shannon, x=pH, group=Buffer)) +
  geom_line(aes(linetype=Buffer)) +
  geom_point(size=1) +
  ylab('Mean Shannnon diversity of replicates') +
  theme_minimal () +
  ylim(0,3))
ggsave(paste0(fig_dir, 'subpanels/Fig_1G_shannonDiv.pdf'),width=4,height=3)
