#######################################
# The following script includes analysis that led to the creation of the synthetic communities
# 
# Rebecca Culver
#######################################

# configure directories, load libraries and base functions
rm(list=ls())
source(paste0(here::here(), "/0-config.R"))
ps<-readRDS(paste0(clean_data_dir,'merge_phyloseq.RDS'))
genus.order<-readRDS(paste0(fig_dir, 'official_colors.RDS'))

#######################################
# In text discussion of isolate library
#######################################

iso.df <- read.csv(paste0(raw_data_dir,'220825_isolateData.csv'))

unique(iso.df$Sample_Type)

print('Location of aquisition')
iso.df %>% group_by(Sample_Type) %>% dplyr::summarise(n())

iso.df.noHumanStool <- iso.df %>% filter(!Sample_Type %in% 'humanStool')
dim(iso.df.noHumanStool)

print('Number of Phyla')
length(unique(iso.df.noHumanStool$Phylum))
print('Number of Families')
length(unique(iso.df.noHumanStool$Family))
print('Number of Species')
length(unique(iso.df.noHumanStool$Species))

# Families in isolate library that were also in the human intestinal sample
df2plot<- ps %>% 
  subset_samples(SampleType %in% c('original_2522_capsule') | Project %in% 'c2m5' & pH == 0 & !SampleType %in% 'regrow') %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>%
  filter(Abundance > 0.001) %>%
  mutate(Inoculum = ifelse(is.na(Inoculum)==TRUE, '2522', Inoculum))
unique(df2plot$Inoculum)

print('Abundance represented in devices by Family')
df2plot %>% group_by(Inoculum) %>% filter(Family %in% iso.df$Family) %>% dplyr::summarise(total=sum(Abundance)) %>% arrange(total)
#Sample 2522 was used to generate SICom


# Families in islate library that were also in invitro communities
ps.tmp <- ps
sample_data(ps.tmp) <- data.frame(sample_data(ps.tmp)) %>% mutate(id= paste0(Inoculum, '_',pH,'_', Buffer))
ps2plot<- ps.tmp %>% 
  subset_samples(SampleType %in% c('DAY0_1') | Project %in% 'c2m5' & !pH==0 & !SampleType %in% 'regrow') %>%
  merge_samples("id")%>%
  transform_sample_counts(function(x) x/sum(x)) 

df2plot <- ps2plot %>%
  psmelt() %>% data.frame() %>%
  filter(Abundance> 0.001)

print('Abundance represented in in vitro communities by Family')
df2plot %>% filter(Family %in% iso.df$Family) %>% group_by (Sample) %>% 
  dplyr::summarise(Abundance=sum(Abundance)) %>% arrange(Abundance) %>% data.frame()
#2522_7.5_100 is the SICom, StoolCom is 197_7_NA 
# The NA at the end of all the sample names in this data frame refers to 0 added buffer


# df of small intestine isolates from human ONLY
iso.siily <- iso.df %>% filter(!Sample_Type %in% c('humanStool','mouseStool','mouseSmallIntestine'))
dim(iso.siily)

# What did the mouse colonization accomplish...
print('MouseStool not initially present')
mouseStool <- iso.df %>% filter(Sample_Type %in% 'mouseStool')
dim(mouseStool)
dim(mouseStool %>% filter(!Species %in% iso.siily$Species))

print('MouseSmallIntestine not initially present')
mouseSi <- iso.df %>% filter(Sample_Type %in% 'mouseSmallIntestine')
dim(mouseSi)
dim(mouseSi %>% filter(!Species %in% iso.siily$Species))

## Comparing humanStool strains to what was acquired from the small intestine

print('HumanStool unique species')
humanStool <- iso.df %>% filter(Sample_Type %in% 'humanStool')
dim(humanStool)
length(unique(humanStool$Species))
length(unique(humanStool$Family))
length(unique(humanStool$Phylum))

print('HumanStool species not found in library')
humanStool <- iso.df %>% filter(Sample_Type %in% 'humanStool')

dim(humanStool %>% filter(!Species %in% iso.siily$Species))
#humanStool %>% filter(!Species %in% iso.siily$Species)
dim(humanStool %>% filter(!Family %in% iso.siily$Family))
humanStool %>% filter(!Family %in% iso.siily$Family)


#######################################
# Figure 4A
#######################################


iso.df.plot <- iso.df %>%
  group_by(Family) %>%
  dplyr::summarise(iso_num = n()) %>% 
  melt(id=c('iso_num')) %>%
  dplyr::rename(Family=value) %>% select(-variable) %>%
  left_join(genus.order, by='Family') %>%
  mutate(Family = factor(Family, levels=genus.order$Family)) %>%
  arrange(Phylum, Family)

p<-ggplot(iso.df.plot, aes(x='', y=iso_num, fill=Family)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=iso.df.plot$colors) +
  coord_polar("y") +
  theme_minimal()
ggsave(paste0(fig_dir, 'subpanels/Fig_4A_isolates.pdf'),p, width=6, height=6)


