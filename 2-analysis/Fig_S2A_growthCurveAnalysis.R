#######################################
# The following script includes analysis for the
#  growth curves of in vitro communities
# BCECF was used to track pH change
# OD600, 440, and 485 were measured
# 
# Rebecca Culver
#######################################

# configure directories, load libraries and base functions
rm(list=ls())
source(paste0(here::here(), "/0-config.R"))


######################################
# The following imports the od data and metadata for each sample
# It also calibrates the 440/485 readings to give the pH
######################################

# Variables
vars <- read.csv(paste0(raw_data_dir,'plateReader/vars_mp1_mp2.csv')) %>%
  dplyr::rename(Media_Plate=Meda_Plate)

# Media plate 1 c2m6
f<-read.csv(paste0(raw_data_dir,'plateReader/201217_c2m6_d10_p3_mp1_BCECF.csv')) 

# Dataframe restructuring
df<-melt(f, id=c("Time","Temp","Read")) 
colnames(df)<-c('Time', 'Temp', "Read", 'Well','OD')
df$Time<-as.matrix(read.table(text = df$Time, sep = ":")) %*% c(60, 1, 1/60) #Convert time to hours
df$Time<-df$Time/60
df <- df %>% 
  left_join(vars %>% filter(Media_Plate %in% 'mp1'), by='Well') %>%
  filter(!Well %in% c('H1','H2','H3')) %>% # Remove wells were buffer was done twice (pH6.5 -> MES, removed MOPS)
  mutate(c2m = 'c2m6')
mp1 <- df


# Media plate 2 c2m6
# Add in variable info about samples
f<-read.csv(paste0(raw_data_dir,'plateReader/210111_c2m6_p3_mp2_BCECF.csv'))

# Dataframe restructuring
df<-melt(f, id=c("Time","Temp","Read")) 
colnames(df)<-c('Time', 'Temp', "Read", 'Well','OD')
df$Time<-as.matrix(read.table(text = df$Time, sep = ":")) %*% c(60, 1, 1/60) #Convert time to hours
df$Time<-df$Time/60
df <- df %>% 
  left_join(vars %>% filter(Media_Plate %in% 'mp2'), by='Well') %>%
  mutate(c2m = 'c2m6') %>%
  filter(!Well %in% c('E1','E2','E3','C1','C2','C3'))  # Remove wells were buffer was done twice (pH 7.5 -> MOPS, remove HEPES; pH 8 -> HEPES, removed)
mp2 <- df



# Add in the 100 mM buffer conditions
f<-read.csv(paste0(raw_data_dir,'plateReader/210116_c2m7_d14_p4_mp1_bcecf.csv'))

# Dataframe restructuring
df<-melt(f, id=c("Time","Temp","Read")) 
colnames(df)<-c('Time', 'Temp', "Read", 'Well','OD')
df$Time<-as.matrix(read.table(text = df$Time, sep = ":")) %*% c(60, 1, 1/60) #Convert time to hours
df$Time<-df$Time/60
df <- df %>% 
  left_join(vars %>% filter(Media_Plate %in% 'mp1'), by='Well') %>%
  mutate(Buffer = ifelse(Buffer == 50, 100, Buffer))%>%
  filter(!Well %in% c('G1','G2','G3')) %>% # Remove wells were buffer was done twice (pH6.5 -> MOPS, removed MES)
  mutate(c2m = 'c2m7')
mp1_100 <- df

f<-read.csv(paste0(raw_data_dir,'plateReader/210128_c2m7_d14_p4_mp2_bcecf.csv'))

# Dataframe restructuring
df<-melt(f, id=c("Time","Temp","Read")) 
colnames(df)<-c('Time', 'Temp', "Read", 'Well','OD')
df$Time<-as.matrix(read.table(text = df$Time, sep = ":")) %*% c(60, 1, 1/60) #Convert time to hours
df$Time<-df$Time/60
df <- df %>% 
  left_join(vars %>% filter(Media_Plate %in% 'mp2'), by='Well') %>%
  mutate(Buffer = ifelse(Buffer == 50, 100, Buffer)) %>%
  filter(!Buffer==0) %>%
  filter(!Well %in% c('E1','E2','E3','C1','C2','C3','A2')) %>% # Remove wells were buffer was done twice (pH 7.5 -> MOPS, remove HEPES; pH 8 -> HEPES, removed)
  mutate(c2m = 'c2m7')
mp2_100 <- df

# Create the final dataframe
df.final <- rbind(mp1, mp2, mp1_100, mp2_100) %>%
  filter(Inoculum %in% c('2522','0')) %>%
  mutate(Read = substr(Read,1,3)) %>%
  select(-Temp) %>%
  group_by(c2m, Read, Well, Media_Plate, pH, Buffer, Inoculum) %>%
  mutate(timePosition = 1:n())

time <- df.final %>% select(Time, timePosition, c2m, Read, Well, Media_Plate, pH, Buffer, Inoculum) %>%
  filter(Read == 600) %>% select(-Read)

df.noTime <- df.final %>%
  select(-Time) %>% pivot_wider(names_from = Read, values_from = OD) 

df.final <- df.noTime %>% 
  left_join(time, by=c('timePosition', 'c2m', 'Well', 'Media_Plate', 'pH', 'Buffer', 'Inoculum')) %>%
  select(-Read) %>%
  dplyr::rename(x600='600',x440='440',x485='485')

# Calibration of 440/485 for pH measurements
calibrate.c2m6 <- df.final %>% 
  filter(timePosition==1, Inoculum %in% c(1,2522)) %>%
  filter(c2m=='c2m6') %>%
  mutate(ratio=x440/x485)

p<-ggplot(calibrate.c2m6, aes(y=pH, x=x440/x485, fill=Time)) +
  geom_point() 

lm(pH ~ ratio, calibrate.c2m6)


calibrate.c2m7 <- df.final %>% 
  filter(timePosition==1, Inoculum %in% c(1,2522)) %>%
  filter(c2m=='c2m7') %>%
  mutate(ratio=x440/x485)

p<-ggplot(calibrate.c2m7, aes(y=pH, x=x440/x485, fill=Time)) +
  geom_point() 

lm(pH ~ ratio, calibrate.c2m7)

# Now add this calibration from the above lm to the dataframe
df.final <- df.final %>%
  mutate(ratio = x440/x485) %>%
  mutate(pH_measured = ifelse(c2m == 'c2m6', -14.1*ratio + 10.4,
                              -42.9*ratio + 10.1))

######################################
# Generate the plot of the pH min/max during growth
######################################

df2plot <- df.final %>% filter(Inoculum==2522, Time<36, !pH %in% 5.5) %>%
  mutate(to_exclude = ifelse(c2m=='c2m7' & Buffer==0, 1,0)) %>%
  mutate(Buffer= as.character(Buffer)) %>%
  mutate(Buffer = factor(Buffer, levels=c('0','50','100'))) %>%
  filter(!to_exclude==1) %>%
  group_by(Well, pH, Buffer) %>%
  dplyr::summarise(`Minimum pH` = min(pH_measured), `Maximum pH`=max(pH_measured)) %>%
  melt(id=c('Well','pH','Buffer')) %>% ungroup()

ggplot(df2plot, aes(y=value, x=Buffer, shape=variable)) +
  geom_point() +
  theme_minimal() +
  scale_shape_manual(values=c(1,15))+
  ylab('Min and max pH during growth') +
  ggtitle('Starting pH of media') +
  facet_wrap(~pH, nrow=1)
ggsave(paste0(fig_dir, 'subpanels/Fig_S2A_growthCurveAnalysis.pdf'), width=7,height=3)
