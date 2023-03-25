#######################################
# The following script includes analysis for bottom-up community assembly
# 
# Rebecca Culver
#######################################

# configure directories, load libraries and base functions
rm(list=ls())
source(paste0(here::here(), "/0-config.R"))
ps<-readRDS(paste0(clean_data_dir,'merge_phyloseq.RDS'))
genus.order<-readRDS(paste0(fig_dir, 'official_colors.RDS'))
species.order<-read.csv(paste0(raw_data_dir, 'species_colors.csv'))

#######################################
# Create a phyloseq object that merges all the SI samples together and reannotates dataframe
#######################################

ps.synMice<- ps %>%
  subset_samples(Project %in% c('MouseExp3') & finalDaySamples %in% 1 & !SampleType %in% c('Cecum','original_2522_capsule')) %>% 
  merge_samples("mouse_location") %>%
  prune_taxa(taxa_sums(.) > 0, .) 
sam_names<-rownames(sample_data(ps.synMice))                            
sample_data(ps.synMice) <- data.frame(ps.synMice@sam_data) %>%  
  mutate(SeqName = sam_names) %>%
  mutate(Mouse = str_split(SeqName, '_', simplify=TRUE)[,1],
         Location = str_split(SeqName, '_', simplify=TRUE)[,2],
         Diet = str_split(SeqName, '_', simplify=TRUE)[,3],
         Community = str_split(SeqName, '_', simplify=TRUE)[,4],
         Community = str_replace_all(Community, ' ','_'),
         mice_loc = paste0(Location,'_', Community)) %>%
  mutate(Diet = factor(Diet, levels=c('Normal','noFiber'))) 

#######################################
# Fig. 4B: Generation of the tree used in the panel that was further edited/colorized in illustrator
#######################################

ps.i<-ps %>%
  subset_samples(Project %in% c('m3g1','m3g2') & SampleType %in% 'Day0_Gavage') %>%
  filter_taxa(function(x) sum(x) > 10, TRUE) %>% tip_glom(h=0.1)
tax_table(ps.i)
new_names <- data.frame(tax_table(ps.i)) %>% 
  rownames_to_column("ASV") %>%
  mutate(Species = ifelse(is.na(Species)==TRUE, ASV, Species)) %>%
  mutate(name=paste0(Genus, ' ',Species))
taxa_names(ps.i) <- new_names$name


pdf(paste0(fig_dir, 'subpanels/Fig_4B_tree_of_inoculated_species.pdf'), width=26, height=6)
plot_tree(ps.i, "treeonly",ladderize='right',label.tips="taxa_names")# + coord_polar(theta="y")
dev.off()

#######################################
# Supplementary Table 4: SI enrichment analysis used to inform generation of SISynCom
#######################################

mouse.ps <- ps %>% subset_samples(Project %in% c('mouse_exp1') & finalDaySamples %in% 1 & !SampleType %in% 'Cecum') %>%
  subset_samples(Diet %in% 'Normal') %>%
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
# Filter for taxa that are detected in 0.1% abundance in at least 3 mice
taxa2keep <- mouse.ps %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% filter(Abundance > 0.001) %>% group_by(OTU) %>%
  dplyr::summarise(n=n()) %>% filter(n>2)

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
  select(OrgName, logFC, adj.P.Val,col) %>% arrange(col, logFC)

condensed_si_microbes <- limmaRes_sig %>%
  filter(col %in% 'Increased in SI') %>%
  mutate(Genus = str_split(OrgName, ' ', simplify=TRUE)[,1],
         Species = str_split(OrgName, ' ', simplify=TRUE)[,2]) %>%
  mutate(species_name = paste0(Genus,' ', Species)) %>%
  select(species_name, logFC, adj.P.Val,col) %>% unique()
write.csv(condensed_si_microbes, paste0(fig_dir,'tables_unannotated/si_enriched_microbes_for_synSI_creation.csv'), row.names=FALSE)


#######################################
# Supplementary Table 5: All isolates that were input into the mice
#######################################

# The first gavage:
# There were no FullSynCom so I just did a spot check to make sure everything checks out"
# also be aware that the isolates were in a random order because I enjoy pain apparently, the "Community" is the isolate it is supposed to be
gavage1 <- ps %>% subset_samples(Project %in% 'm3g1' & Mouse %in% 'Isolate')

# Okay looking at the second gavage which contains all communities to construct this suppolemental table
# the "Community" is the isolate it is supposed to be
gavage2 <- ps %>% subset_samples(Project %in% 'm3g2' & !SampleType %in% 'Day0_Gavage' & !Community %in% c('EMPTY','')) %>%
  filter_taxa(function(x) sum(x) > 1000, TRUE) %>% # because these are pure isolates %>%
  transform_sample_counts(function(x) x/sum(x))

# Check for contamination
table<-gavage2 %>% psmelt() %>% filter(Abundance > 0.1) %>% select(SeqName, Community, Abundance, Genus, Species) %>% arrange(SeqName)

# We can see that a few contaminants:
# m3g2_p1_F12_well_A8_Isolate_S72 - Supposed to be Bifidobacterium but it's actually a L rhamnosus culture
# m3g2_p1_F5_well_A1_Isolate_S65 - Phascolarctobacterium contains a Lactococcus species
# m3g2_p1_G12_well_B8_Isolate_S84 - Slackia culture contains a Lactococcus species
# m3g2_p1_H4_well_B12_Isolate_S88 - Supposed to be a Streptococcus species. This strian did not grow well in the overnight culture, so my guess is that this is cross-contamination from other wells
# m3g2_p1_H8_well_C4_Isolate_S92 - Similar to the aboce Streptococcus species, Colliinsella does not grow well, so we're liekly seeing a lot of cross-contamination

# Everything else looks good - output to a csv that we'll go on to further annotate in excel
write.csv(table, paste0(fig_dir, 'tables_unannotated/isolates_for_gavage.csv') )


#######################################
# Fig. 4C: Mouse experiment schematic overview made in adobe illustrator
#######################################

#######################################
# Fig. 4D: Abundance of different taxa in the synthetic community cross-colonization
# Geom-tile of abundances highlights priority effects
#######################################

ps2plot <- ps.synMice %>% 
  subset_samples(Community %in% c('synIC','synSC_into_synIC-mouse','synSC','synIC_into_synSC-mouse') &
                 Diet %in% 'Normal') %>%
  merge_samples("mice_loc")

taxa2keep<-psmelt(ps2plot %>% transform_sample_counts(function(x) x/sum(x))) %>%
  filter(Abundance > 0.001) %>% group_by(Community, Location, Diet, OTU) %>% 
  dplyr::summarise(num_OTU = n()) %>%
  filter(num_OTU > 1) # Make it for at least two samples

ps2plot <- prune_taxa(taxa2keep$OTU, ps2plot)

df <- psmelt(ps2plot %>% transform_sample_counts(function(x) x/sum(x))) %>%  
  mutate(Abundance = ifelse(Abundance < 0.001, 0.001, Abundance)) %>%
  mutate(Abundance = log10(Abundance)) %>%
  mutate(Location = str_split(Sample, '_',simplify=TRUE)[,1])%>%
  mutate(Sample = factor(Sample, levels=c('Small Intestine_synIC','Small Intestine_synSC_into_synIC-mouse','Small Intestine_synSC','Small Intestine_synIC_into_synSC-mouse',
                                          'Stool_synIC','Stool_synSC_into_synIC-mouse','Stool_synSC','Stool_synIC_into_synSC-mouse'))) %>%
  mutate(Species = ifelse(is.na(Species)==TRUE, 'sp', Species)) %>%
  mutate(Species_name = paste0(Genus, ' ', Species)) %>%
  mutate(Species_name2 = paste0(substr(Genus,1,1),'. ',Species))


ggplot(df, aes(x=Sample, y=Species_name, fill=Abundance)) +
  geom_tile() +
  scale_fill_gradient(low='white',high='blue') +
  theme_minimal() 
ggsave(paste0(fig_dir,'subpanels/Fig_4D_abundanceCross.pdf'), width=6, height=4)

#######################################
# Fig. 4E: Taxa ordination
#######################################

mouse.ps <-ps.synMice %>% 
  subset_samples(Community %in% c('synIC','synSC_into_synIC-mouse','synSC','synIC_into_synSC-mouse','synAC-mouse') &
                   Diet %in% 'Normal') 

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
  mutate(Diet = ifelse(Diet %in% 'noFiber','MD','SD')) %>%
  mutate(Diet = factor(Diet, levels=c('SD','MD'))) 

taxa2keep <- mouse.ps %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% filter(Abundance > 0.001)

mouse.ord.ps <- prune_taxa(taxa2keep$OTU, mouse.ps)  %>% transform_sample_counts(function(x) x/sum(x))

pcoa.bray <- ordinate(mouse.ord.ps,  method = "PCoA", distance = "bray")

var_exp <- get_evals(pcoa.bray)$variance_exp
scores <- get_scores(pcoa.bray, sample_data(mouse.ord.ps))

plot_ordination(mouse.ord.ps, pcoa.bray, color="Community",shape="Community") +
    geom_point(size=1)+
    scale_color_manual(values=c('black','blue','plum','red','brown4'))+
    scale_shape_manual(values=c(1,1,1,1,1))+
    theme_minimal()+
    coord_fixed(sqrt(var_exp[2] / var_exp[1])) +
    facet_wrap(~Location, nrow=1)
ggsave(paste0(fig_dir,'subpanels/Fig_4E_bray_ordination.pdf'), width=8, height = 3)


#######################################
# Fig. 4F,G,H: Taxa maintain spatial niche along intestinal tract during bottom-up community colonization
#######################################
#'synSC_into_synIC-mouse' and 'synIC_into_synSC-mouse' and 'synFullCom'(synAC-mouse) were all generated using the same script
#'The community parameter was swapped out based on which communityy was being analyzed

taxa2keep<-psmelt(ps.synMice %>% 
                    subset_samples(Community %in% c('synAC-mouse') & Diet %in% 'Normal') %>%
                    transform_sample_counts(function(x) x/sum(x))) %>%
  filter(Abundance > 0.001) %>% group_by(OTU) %>% dplyr::summarise(count=n()) %>% filter(count>2)

ps_da <- ps.synMice %>% subset_samples(Community %in% c('synAC-mouse') & Diet %in% 'Normal')
sample_data(ps_da) <- data.frame(sample_data(ps_da)) %>% mutate(mouse_subj = paste0(Mouse,'_',MouseType,'_',Community))

ps_da <- prune_taxa(taxa2keep$OTU, ps_da)


## Estimate size factors
dds <- phyloseq_to_deseq2(ps_da, design = ~ 1)
dds <- estimateSizeFactors(dds, type = "poscount")
size_factors <- sizeFactors(dds)
summary(size_factors)

# Since we do not want to artificially inflate the counts in voom transformation
# we normalize the size factors so the the smallest is equal to 1.
size_factors <- size_factors/min(size_factors)
summary(size_factors)

# Filter sequences
minSubjPrev <- 2 #since most mice were in cages of 4, let's make it 2

# Filter out rare taxa not present in at minSubjPrev
(pstips <- filter_subject_prevalence(ps_da, thresh = minSubjPrev))

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

require("ggrepel")
limmaRes.df <- limmaRes %>%
  left_join(species.order, by=c('Family','Genus','Species')) %>%
  #mutate(colors = ifelse(adj.P.Val > alpha, 'black', colors)) %>%
  #mutate(colors = ifelse(abs(logFC) < 0.75, 'black', colors)) %>%
  mutate(Species_name = ifelse(adj.P.Val > alpha, '', Species_name))

ggplot(limmaRes.df, aes(y=-log10(adj.P.Val), x=logFC, label=Species_name)) +
  geom_point(color=limmaRes.df$colors, size=3)+
  geom_hline(yintercept=1.3, linetype="dashed")+
  geom_vline(xintercept=0.75, linetype="dashed")+
  geom_vline(xintercept=-0.75, linetype="dashed")+
  geom_text_repel()+
  theme_minimal()

ggsave(paste0(fig_dir, 'subpanels/Fig_4H_diffEnriched_FullSynCom.pdf'),width=4,height=5)


#######################################
# Text: Phascolarctobacterium abundance
#######################################

taxa2keep<-psmelt(ps.synMice %>% 
                    subset_samples(Community %in% c('synAC-mouse') & Diet %in% 'Normal') %>%
                    transform_sample_counts(function(x) x/sum(x))) %>%
  filter(Abundance > 0.001) %>% group_by(OTU) %>% dplyr::summarise(count=n()) %>% filter(count>2)

phas.df <- ps.synMice %>% 
  subset_samples(Community %in% c('synAC-mouse') & Diet %in% 'Normal') %>% 
  prune_taxa(taxa2keep$OTU, .) %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>%
  filter(Genus %in% 'Phascolarctobacterium') %>%
  group_by(Sample) %>% dplyr::summarise(Abund=sum(Abundance)) %>% data.frame() %>%
  arrange(Abund)
phas.df

#######################################
# Fig. 4I,J: Ordinations of functional differences along the intestinal tract
#######################################

# This will take a minute
pfam.df<-readRDS(paste0(raw_data_dir,'pfam_output/df_final.RDS'))
pfam.filt <- pfam.df %>% filter(!Community %in% c('Complete','Small_Intestine','Stool',''), 
                                !Location %in% 'Cecum', Diet %in% 'Normal') %>%
  filter(!BioHub_Name %in% 'MouseExp3_p2_B12_Distal_m8Normal') # This mouse received on one SI sample and it's a major outlier in the dataset

# Note that due to poor sequencing depth, the proximal, mid, and distal intestine of a single mouse
# is summed together to give one small intestine sample per mouse

# We will be normalizing each sample to the total number of Pfams detected

#Bonus info: sequencing depth
seq.depth <- pfam.filt %>% select(Mouse, Community, Location, hostRemoved) %>%
  unique() %>% group_by(Mouse, Community, Location) %>%
  dplyr::summarise(hostRemoved = sum(hostRemoved))

pfam.filt <- pfam.df %>% filter(!Community %in% c('Complete','Small_Intestine','Stool',''))
                                
total.pfams <- pfam.filt %>%
  select(Mouse, Community, Location, Diet, pfam) %>% 
  group_by(Mouse, Community, Location, Diet) %>%
  dplyr::summarise(total_pfams = n())

# Now we want to create a counts table for all Pfam samples
pivot.pfam <- pfam.filt %>% group_by(Mouse, Community, Location, Diet, pfam) %>%
  dplyr::summarise(pfam_count = n()) %>% ungroup() %>%
  #left_join(seq.depth, by=c('Mouse', 'Community', 'Location')) %>%
  left_join(total.pfams, by=c('Mouse','Community','Location','Diet')) %>% 
  #mutate(pfam_norm_count = (pfam_count/hostRemoved)*1000000) %>%
  mutate(pfam_norm_count = (pfam_count/total_pfams)*100) %>%
  select(Mouse, Community, Location, pfam, pfam_norm_count) %>%
  pivot_wider(values_from='pfam_norm_count', names_from='pfam', values_fill=0) %>% data.frame()
pivot.pfam.matrix <- pivot.pfam[,-c(1,2,3)]

# PCA
res.pca <- prcomp(pivot.pfam.matrix)
# Eigenvalues
eig.val <- get_eigenvalue(res.pca)
# Results for Variables
res.var <- get_pca_var(res.pca)
# Results for individuals
res.ind <- get_pca_ind(res.pca)

#res.X$coord          # Coordinates
#res.X$contrib        # Contributions to the PCs
#res.X$cos2           # Quality of representation 

pca.ind.df <- as.data.frame(res.ind$coord) %>% select(Dim.1, Dim.2) %>% 
  cbind(as.data.frame(pivot.pfam[,c(1,2,3)]))

ggplot(pca.ind.df, aes(y=Dim.2, x=Dim.1, color=Community, shape=Community)) +
  geom_point(size=1) + 
  scale_color_manual(values=c('blue','brown4','black','red','plum'))+
  scale_shape_manual(values=c(1,1,1,1,1))+
  xlab(paste0('Principal component 1\n(',round(eig.val$variance.percent[1],0),'%)')) +
  ylab(paste0('Principal component 2\n(',round(eig.val$variance.percent[2],0),'%)')) +
  # typically would represent an ordination as a proportion of the variance
  # however, the PCA2 is quite small here
  #coord_fixed(sqrt(eig.val$variance.percent[2] /eig.val$variance.percent[1])) +
  # so instead,
  coord_fixed(1.5/1)+
  facet_wrap(~Location, nrow=2)+
  theme_minimal()

ggsave(paste0(fig_dir, 'subpanels/Fig_4I_pfam_pca.pdf'),width=5,height=3)

var.contributions <- as.data.frame(res.var$contrib) %>% select(Dim.1,Dim.2)%>% 
  rownames_to_column("pfam") %>% 
  dplyr::rename(contribution1=Dim.1, contribution2 = Dim.2)
pca.var.df <- as.data.frame(res.var$coord) %>% select(Dim.1, Dim.2) %>% 
  rownames_to_column("pfam") %>% 
  left_join(var.contributions, by='pfam') %>% arrange(contribution1, decreasing=TRUE)

top.var.contributors <- pca.var.df[1:10,]

ggplot(top.var.contributors, aes(x=Dim.1, y=Dim.2, label=pfam)) +
  geom_point(fill='white', color='white') +
  geom_segment(aes(xend=0, yend=0), arrow = arrow(length = unit(0.1,"cm"), ends="first")) +
  theme_minimal() +
  geom_text_repel() +
  coord_fixed(1.5/1)

ggsave(paste0(fig_dir, 'subpanels/Fig_4J_pfamVars_pca.pdf'),width=8,height=6)

## Top pfam definitions
# CarbopepD_reg_2 = mediates substrate-specific transport
# TonB_dep_Rec = TonB dependent receptor; mediates substrate-specific transport
# ABC_tran = ABC transporters for a large family of proteins responsible for translocation of a variety of compounds across biological membranes.
# BPD_transp_1 = transmembrane protein of unknown function
# OEP = Outer membrane efflux protein that transports substrates
# Plug = blocks channel of a receptor until bound by ligand
# SusD.like_3 = SusD a secreted poly-saccharide-binding protein. Probably mediates xyloglucan-binding prior to xyloglucan transport in the periplasm for degradation
# SusD_RagB = Structural characterization of RagB shows substantial similarity with Bacteroides thetaiotaomicron SusD (i.e alpha-helices and TPR regions). Based on this structural similarity, functional studies suggest that, 
# RagB binding of glycans occurs at pockets on the molecular surface that are distinct from those of SusD 

# Reg_prop = two-component regulator protein; form a beta-propeller; unclear function
# HTH_18 = domain associated with AraC, which encodes a positive regulator protein for L-arabinose utilization

# 2 of these proteins is a regulator (1 for L arabinose utilization)
# 8 of these proteins are known (or putatively predicted) to be iunvolved with substrate import/export;
#         1 of these is a plug, two of these are SusD-like and are secreted


# How many unique pfams are in each of these samples? Bottom-up synthetic communities on normal diet
m3.normal <- pfam.df %>% filter(Project %in% 'MouseExp3', Diet %in% 'Normal', Community %in% 'syn_SI_syn_Stool_syn_All')
#. Sequencing depth analysis
m3.normal %>% select(BioHub_Name, Community, Location, hostRemoved) %>% unique() %>%
  group_by(Location) %>% dplyr::summarise(sum(hostRemoved))

print('Total')
m3.normal %>% select(pfam) %>% unique() %>% dplyr::summarise(unique_pfams=n())
print('Small intestine')
m3.si <- m3.normal %>% filter(Location %in% 'Small Intestine') %>% select(pfam) %>% unique() 
num<-m3.normal %>% filter(Location %in% 'Small Intestine') %>% select(pfam) %>% unique() %>% dplyr::summarise(unique_pfams=n())
num$unique_pfams

print('Stool')
m3.stool<- m3.normal %>% filter(Location %in% 'Stool') %>% select(pfam) %>% unique() 
num<-m3.normal %>% filter(Location %in% 'Stool') %>% select(pfam) %>% unique() %>% dplyr::summarise(unique_pfams=n())
num$unique_pfams

# How many were only found in the small intestine or stool?
print('Stool only')
dim(m3.stool %>% filter(!pfam %in% m3.si$pfam))
print('SI only')
dim(m3.si %>% filter(!pfam %in% m3.stool$pfam))

# Unique to SI and found in all 4 mice
all.mice.si<-m3.normal %>% filter(Location %in% 'Small Intestine', !pfam %in% m3.stool$pfam) %>% 
  group_by(Mouse, pfam) %>% dplyr::summarise(pfam_count = n()) %>% 
  mutate(pfam_count = ifelse(pfam_count > 0,1,0)) %>% group_by(pfam) %>% 
  dplyr::summarise(mouse_count=n()) %>%
  filter(mouse_count > 3)

write.csv(all.mice.si, paste0(fig_dir, 'tables_unannotated/unique_pfams_in_si.csv'))

### Bonus: Creating pfam dataframe. Please note that this script takes a while to run. Therefore,
# after it was run once, it was simply saved and reloaded as the df_final.RDS

vars <- read.csv(paste0(raw_data_dir,"metagenomics_organized/variables_metaData.csv")) %>% select(-X)
files <- list.files(path=paste0(raw_data_dir, 'pfam_output'),pattern="*resolved", full.names=TRUE)
genome.df <- lapply(files, function(x) read.delim(x, skip=2, sep=' ', header=FALSE))
genome.df <- rbind.fill(genome.df) %>%
  dplyr::rename(query=V1, pfam=V2, score=V3, boundaries=V4, resolved=V5, cond_evalue=V6, indp_evalue=V7)

df.final <- genome.df %>%
  mutate(BioHub_Name = str_split(query, '.scaffolds', simplify=TRUE)[,1]) %>%
  mutate(BioHub_Name = str_split(BioHub_Name, '_S',simplify=TRUE)[,1]) %>%
  mutate(scaffold = str_split(query, '.scaffolds_', simplify=TRUE)[,2]) %>%
  left_join(vars, by='BioHub_Name') %>%
  # remove contaminated stool samples
  filter(!BioHub_Name %in% c('mousexp3_gavage1_D2','mousexp3_gavage1_D3','mousexp3_gavage1_B1')) %>%
  filter(hostRemoved > 1000000) %>%# filter low read samples
  mutate(Location = ifelse(!SampleType %in% c('Cecum','Stool'), 'Small Intestine', SampleType )) 

saveRDS(df.final, paste0(raw_data_dir,'pfam_output/df_final.RDS'))
