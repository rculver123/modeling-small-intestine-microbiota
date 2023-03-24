#######################################
# The following script includes analysis for the
#  growth curves for Enterococcus isolates
# 
# Rebecca Culver
#######################################

# configure directories, load libraries and base functions
rm(list=ls())
source(paste0(here::here(), "/0-config.R"))


########################### 

# Variables
vars <- read.csv(paste0(raw_data_dir,'entero_CarbonSource_Screen/entero_screen_vars.csv')) %>%
  mutate(id=paste0(plateNum,Well)) 

f<-list.files(paste0(raw_data_dir,'entero_CarbonSource_Screen/plateReader'), full.names=TRUE) 

all_files_merged <- data.frame()
for (item in f) {
  name <- gsub("^.*/", "", item)
  name_stripped <- gsub(".csv","",name)
  plate_num <- str_split(name_stripped, '_', simplify=TRUE)[,3]
  tmp.df <- read.csv(item) %>%
    melt(id=c("Time","Temp")) %>%
    mutate(plateNum = plate_num)
  all_files_merged <- rbind(all_files_merged, tmp.df)
}

# Convert time to hours
df<-all_files_merged
df$Time<-as.matrix(read.table(text = df$Time, sep = ":")) %*% c(60, 1, 1/60) #Convert time to hours
df$Time<-df$Time/60

# Create final dataframe
df.final <- df %>%
  dplyr::rename(Well=variable, od600=value) %>%
  left_join(vars, by=c('Well','plateNum')) %>%
  filter(Time < 48) # things began evaporating after this

# Calculate the area under the curve (auc) for each of the samples
grouped.list<-split(df.final, df.final$id)
# Helper function to run Growthcurver's SummarizeGrowth
get_auc <- function(df) {
  summary <- SummarizeGrowth(as.vector(df$Time), as.vector(df$od600))
  auc <- summary$vals$auc_e
  return(auc)  
}
auc_e <- lapply(grouped.list, function(x) get_auc(x)) 
auc.df <- t(data.frame(auc_e)) %>% as.data.frame() %>%
  dplyr::rename(auc=V1) %>% rownames_to_column('id') %>%
  right_join(vars, by='id')

# Figure out if there are any crazy outliers (potentially that didn't grow)
hist(auc.df$auc, breaks=25)
df.sort <- auc.df %>% arrange(auc)
top_outliers <- df.sort[1:3,]
# Plots of outliers
df.outliers <- df.final %>%
  filter(id %in% top_outliers$id)
(p<-ggplot(df.outliers, aes(y=(od600), x=Time, color=Strain )) +
    geom_point(size=0.5) +
    theme_minimal() +
    facet_wrap(~Carbon))
# p4G10 appears to be an outlier that did not grow


# Create final data frame to plot
df.auc.plot <- auc.df %>% filter(!id %in% 'p4G10') %>% # filter out this well that didn't grow
  mutate(Type = ifelse(plateNum %in% c('p1','p2'), 'Simple','Complex')) %>%
  mutate(Type = factor(Type, levels=c('Simple','Complex'))) %>%
  group_by(Strain, Carbon, Type) %>%
  dplyr::summarise(mean_auc = mean(auc), sd_auc=sd(auc)) %>% ungroup()

# Visualize std deviations- everything checks out
histogram(df.auc.plot$sd_auc)

# Normalize complex carbohydrates to the glucose only condition
glucose.only <- df.auc.plot %>% filter(Carbon %in% 'Glucose') %>%
  dplyr::rename(glu_auc=mean_auc) %>% select(Strain, glu_auc)
df.auc.norm <- df.auc.plot %>%
  left_join(glucose.only, by='Strain') %>%
  mutate(mean_auc = ifelse(Type %in% 'Simple', mean_auc, mean_auc/glu_auc)) %>%
  mutate(Strain = factor(Strain, levels=c("Enterococcus_DualCom_C","Enterococcus_DualCom_A","Enterococcus_DualCom_B",
                                          "Enterococcus_StoolCom_A","Enterococcus_StoolCom_B","Enterococcus_SynStoolCom")))

# Plot simple carbon sources
ggplot(df.auc.norm %>% filter(!Type %in% 'Simple'), aes(x=Carbon, y=Strain, fill = mean_auc)) +
  geom_tile() +
  theme_minimal() +
  scale_fill_gradient(low='white',high='blue') +
  #scale_fill_gradient2(low='white',mid='#6DD5FA',high='blue', midpoint=12.5) +
  theme(axis.text.x = element_text(angle =90)) +
  facet_wrap(~Type, scale='free_x')
ggsave(paste0(fig_dir, 'subpanels/Fig_5X_complex_normalized.pdf'), width=4, height=4)


# Example plot

df2plot <- df.final %>%
  filter(Strain %in% c('Enterococcus_DualCom_C','Enterococcus_StoolCom_A'),
         Carbon %in% c('Sunfiber + Glucose','Inulin + Glucose','PrecticX + Glucose'))

(p<-ggplot(df2plot, aes(y=(od600), x=Time, color=Strain )) +
  geom_point(size=0.5) +
  scale_color_manual(values = c('purple','red'))+
  theme_minimal() +
  facet_wrap(~Carbon) +
  ylim(0,1.4))



######################################
# Old code that generates plots for all samples
######################################



# Read all of the media plates and then merge into one dataframe
f<-list.files(paste0(raw_data_dir,'entero_CarbonSource_Screen/plateReader'), full.names=TRUE) 

all_files_merged <- data.frame()
for (item in f) {
  name <- gsub("^.*/", "", item)
  name_stripped <- gsub(".csv","",name)
  plate_num <- str_split(name_stripped, '_', simplify=TRUE)[,3]
  tmp.df <- read.csv(item) %>%
    melt(id=c("Time","Temp")) %>%
    mutate(plateNum = plate_num)
  all_files_merged <- rbind(all_files_merged, tmp.df)
}

df <- all_files_merged %>% dplyr::rename(Well=variable) %>%
  left_join(vars, by=c('Well','plateNum'))

# Fix time to be in hours
df$Time<-as.matrix(read.table(text = df$Time, sep = ":")) %*% c(60, 1, 1/60) #Convert time to hours
df$Time<-df$Time/60

df.final <- df %>%
  mutate(Type = ifelse(plateNum %in% c('p1','p2'), 'Simple','Complex')) %>%
  mutate(Type = factor(Type, levels=c('Simple','Complex')))


######################################
# Generate plot
######################################

strain_colors <- c('red','red','red','red',
                        'orangered1','orangered1','orangered1','orangered1',
                        'orange','orange','orange','orange',
                        'darkorchid','darkorchid','darkorchid','darkorchid',
                        'turquoise','turquoise','turquoise','turquoise',
                        'black','black','black','black')
df2plot <- df.final %>%
  #filter(Time < 55) %>%# Plate wells started to evaporate after this point
  mutate(Strain_Rep = paste0(Strain,'_',Replicate)) 

p<-ggplot(df2plot, aes(y=(value), x=Time, color=Strain_Rep)) +
  geom_point(size=0.5) +
  scale_color_manual(values = strain_colors)+
  theme_minimal() +
  facet_wrap(~Type+Carbon) +
  ylim(0,1.4)
ggsave(paste0(fig_dir, 'subpanels/Fig_X_enteroGrowth.pdf'), width=14,height=14)
                        

                        
                        
                        
