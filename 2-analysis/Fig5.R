#######################################
# The following script includes analysis for bottom-up community assembly
# folloowed by a diet switch
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
# Fig. 5A: Mouse schematic generated in illustrator
#######################################

#######################################
# Fig. 5B: Mouse schematic generated in illustrator
#######################################

ps.dietSwitch <- ps.synMice %>% subset_samples(Community %in% c('synAC-mouse','synSC_into_synIC-mouse'))

# Calculate taxa to keep
taxa2keep<-psmelt(ps.dietSwitch %>% transform_sample_counts(function(x) x/sum(x))) %>%
  filter(Abundance > 0.001) %>% group_by(Community, Location, Diet, OTU) %>% 
  dplyr::summarise(num_OTU = n()) %>%
  filter(num_OTU > 2) # Make it for at least two samples

df2plot <- psmelt(prune_taxa(taxa2keep$OTU, ps.dietSwitch) %>%
                    transform_sample_counts(function(x) x/sum(x))) %>%
  left_join(species.order, by=c("Family",'Genus','Species')) %>%
  select(SeqName, Abundance, Species_name, Location, Community, Diet, colors) %>%
  mutate(Species_name = factor(Species_name, levels=species.order$Species_name),
         Community = ifelse(Community %in% 'synAC-mouse', 'SynFullCom','SynStoolCom into SynSICom-colonized mice'))
colors2plot <- species.order %>% filter(Species_name %in% df2plot$Species_name)

(b<-ggplot(df2plot , aes(x=SeqName, y=Abundance, fill=Species_name)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=colors2plot$colors) +
  theme_minimal() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  facet_wrap(~Location+Community+Diet, scale='free_x',nrow=2))

ggsave(paste0(fig_dir,'subpanels/Fig_5B_barplotsDietSwitch.pdf'), b,width=7, height=4)

#######################################
# Fig. 5C: Taxa ordination
#######################################

ps.dietSwitch <- ps.synMice %>% subset_samples(Community %in% c('synAC-mouse','synSC_into_synIC-mouse'))

Sample<-rownames(ps.dietSwitch@sam_data)

sample_data(ps.dietSwitch) <- data.frame(ps.dietSwitch@sam_data) %>% 
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

# Calculate taxa to keep
taxa2keep<-psmelt(ps.dietSwitch %>% transform_sample_counts(function(x) x/sum(x))) %>%
  filter(Abundance > 0.001) %>% group_by(Community, Location, Diet, OTU) %>% 
  dplyr::summarise(num_OTU = n()) %>%
  filter(num_OTU > 2) # Make it for at least two samples

mouse.ord.ps <- prune_taxa(taxa2keep$OTU, ps.dietSwitch)  %>% transform_sample_counts(function(x) x/sum(x))

pcoa.bray <- ordinate(mouse.ord.ps,  method = "PCoA", distance = "bray")

var_exp <- get_evals(pcoa.bray)$variance_exp
scores <- get_scores(pcoa.bray, sample_data(mouse.ord.ps))

plot_ordination(mouse.ord.ps, pcoa.bray, color="Community",shape="Diet") +
    geom_point(size=1)+
    scale_color_manual(values=c('black','brown4'))+
    scale_shape_manual(values=c(1,19))+
    theme_minimal()+
    coord_fixed(sqrt(var_exp[2] / var_exp[1])) +
    facet_wrap(~Location, nrow=2)
ggsave(paste0(fig_dir,'subpanels/Fig_5C_bray_ordination.pdf'), width=5, height = 3)


#######################################
# Fig. 5D: CAZymes
#######################################

# Import: metadata of GTDBTK annotation of isolate genomes
vars <- read.table(paste0(raw_data_dir,"drep_siily/gtdbtk.bac120.summary.tsv"), sep="\t", header=TRUE) %>%
  mutate(SeqName = str_split(user_genome, '[.]filtered', simplify=TRUE)[,1]) %>%
  select(user_genome, SeqName, classification) %>%
  mutate(Species = str_split(classification,';s__', simplify=TRUE)[,2]) %>%
  mutate(Species = ifelse(SeqName %in% c('MouseExp3_p3_A3_Stool_m17noFiber_S127.12','m3g2_p1_F5_well_A1_Isolate_S19'), 'Collinsella sp', Species)) %>%
  mutate(Genus = str_split(Species, ' ', simplify=TRUE)[,1]) 
colnames(vars)
# Import: CAZyme annotation
files <- list.files(path=paste0(raw_data_dir,'drep_siily/cazyme'), pattern="*hmmer*", full.names=TRUE)
file <- rbind.fill(lapply(files, function(x) read.delim(x, header=FALSE, sep='\t') %>% mutate(file_name=x)))
headers<-c('hmm_cazyme','hmm_length','query_gene','query_length','evalue','V6','V7','V8','V9','perc_similraity','file_name')
colnames(file) <- headers
df.final<-file %>%
  mutate(short_file_name = str_split(file_name, 'cazyme/', simplify=TRUE)[,2]) %>%
  mutate(SeqName = str_split(short_file_name,'[.]hmmer', simplify=TRUE)[,1]) %>%
  mutate(SeqName = str_split(SeqName, '[.]fasta', simplify=TRUE)[,1]) %>%
  left_join(vars, by="SeqName") %>%
  group_by(SeqName, Species) %>%
  dplyr::summarise(num_cazymes=n()) %>% 
  data.frame() %>%
  arrange(num_cazymes)

# Determine if there are multiple species 
df.final %>% group_by(Species) %>%
  dplyr::summarise(num_sp = n()) %>% filter(num_sp >1)

ggplot(df.final, aes(y=reorder(SeqName,num_cazymes),x=num_cazymes)) +
  geom_bar(stat='identity') +
  scale_y_discrete(labels = df.final$Species) +
  theme_minimal()
ggsave(paste0(fig_dir,'subpanels/Fig_5D_cazymes.pdf'), width=6,height=3)


#######################################
# Fig. 5D,E: Ordinations of functional differences along the intestinal tract
#######################################

# This will take a minute
pfam.df<-readRDS(paste0(raw_data_dir,'pfam_output/df_final.RDS'))
pfam.filt <- pfam.df %>% filter(Community %in% c('syn_SI_syn_Stool', 'syn_SI_syn_Stool_syn_All'), 
                                !Location %in% 'Cecum') %>%
  filter(!BioHub_Name %in% 'MouseExp3_p2_B12_Distal_m8Normal') # This mouse received on one SI sample and it's a major outlier in the dataset

# Note that due to poor sequencing depth, the proximal, mid, and distal intestine of a single mouse
# is summed together to give one small intestine sample per mouse

# We will be normalizing each sample to the total number of Pfams detected

# Bonus info: sequencing depth
seq.depth <- pfam.filt %>% select(Mouse, Community, Location, Diet, hostRemoved) %>%
  unique() %>% group_by(Mouse, Community, Location) %>%
  dplyr::summarise(hostRemoved = sum(hostRemoved))


total.pfams <- pfam.filt %>%
  select(Mouse, Community, Location, Diet, pfam) %>% 
  group_by(Mouse, Community, Location, Diet) %>%
  dplyr::summarise(total_pfams = n())

# Now we want to create a counts table for all Pfam samples
pivot.pfam <- pfam.filt %>% group_by(Mouse, Community, Location, Diet, pfam) %>%
  dplyr::summarise(pfam_count = n()) %>% ungroup() %>%
  left_join(total.pfams, by=c('Mouse','Community','Location','Diet')) %>% 
  mutate(pfam_norm_count = (pfam_count/total_pfams)*100) %>%
  select(Mouse, Community, Location, Diet, pfam, pfam_norm_count) %>%
  pivot_wider(values_from='pfam_norm_count', names_from='pfam', values_fill=0) %>% data.frame()
pivot.pfam.matrix <- pivot.pfam[,-c(1,2,3,4)]

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
  cbind(as.data.frame(pivot.pfam[,c(1,2,3,4)]))

ggplot(pca.ind.df, aes(y=Dim.2, x=Dim.1, color=Community, shape=Diet)) +
  geom_point(size=1) + 
  scale_color_manual(values=c('brown4','black'))+
  scale_shape_manual(values=c(19,1))+
  xlab(paste0('Principal component 1\n(',round(eig.val$variance.percent[1],0),'%)')) +
  ylab(paste0('Principal component 2\n(',round(eig.val$variance.percent[2],0),'%)')) +
  # typically would represent an ordination as a proportion of the variance
  # however, the PCA2 is quite small here
  coord_fixed(sqrt(eig.val$variance.percent[2] /eig.val$variance.percent[1])) +
  # so instead,
  #coord_fixed(1.5/1)+
  facet_wrap(~Location, nrow=2)+
  theme_minimal()

ggsave(paste0(fig_dir, 'subpanels/Fig_5D_pfam_pca.pdf'),width=5,height=3)

# Determine top contributors
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

ggsave(paste0(fig_dir, 'subpanels/Fig_SX_pfamVars_pca_Diet_bothLocations.pdf'),width=8,height=6)

# Now we want to take a closer look at just the small intestine
pivot.pfam <- pfam.filt %>% 
  filter(Location %in% 'Small Intestine') %>%
  group_by(Mouse, Community, Location, Diet, pfam) %>%
  dplyr::summarise(pfam_count = n()) %>% ungroup() %>%
  left_join(total.pfams, by=c('Mouse','Community','Location','Diet')) %>% 
  mutate(pfam_norm_count = (pfam_count/total_pfams)*100) %>%
  select(Mouse, Community, Location, Diet, pfam, pfam_norm_count) %>%
  pivot_wider(values_from='pfam_norm_count', names_from='pfam', values_fill=0) %>% data.frame()
pivot.pfam.matrix <- pivot.pfam[,-c(1,2,3,4)]

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
  cbind(as.data.frame(pivot.pfam[,c(1,2,3,4)]))

ggplot(pca.ind.df, aes(y=Dim.2, x=Dim.1, color=Community, shape=Diet)) +
  geom_point(size=1) + 
  scale_color_manual(values=c('brown4','black'))+
  scale_shape_manual(values=c(19,1))+
  xlab(paste0('Principal component 1\n(',round(eig.val$variance.percent[1],0),'%)')) +
  ylab(paste0('Principal component 2\n(',round(eig.val$variance.percent[2],0),'%)')) +
  # typically would represent an ordination as a proportion of the variance
  # however, the PCA2 is quite small here
  coord_fixed(sqrt(eig.val$variance.percent[2] /eig.val$variance.percent[1])) +
  # so instead,
  #coord_fixed(1.5/1)+
  theme_minimal()

ggsave(paste0(fig_dir, 'subpanels/Fig_5E_pfam_pca_siOnly.pdf'),width=5,height=3)

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
  coord_fixed(sqrt(eig.val$variance.percent[2] /eig.val$variance.percent[1])) 

ggsave(paste0(fig_dir, 'subpanels/Fig_5F_pfamVars_pca_Diet_siOnly.pdf'),width=8,height=6)

# ACR_tran = membrane protein complex; some involved with drug resistance; 
# https://rdcu.be/c7kxk
