#Quantification analysis for RNAscope measurements 
library(dplyr)
library(ggplot2)
#library(ggpubr)
library(here)
library(ggsci)

# make colour palette
mypal <- pal_npg("nrc", alpha = 0.7)(10)
mypal2 <-pal_tron("legacy", alpha = 0.7)(7)
mypal3 <- pal_lancet("lanonc", alpha = 0.7)(9)
mypal4 <- pal_simpsons(palette = c("springfield"), alpha = 0.7)(16)
mypal5 <- pal_rickandmorty(palette = c("schwifty"), alpha = 0.7)(6)
mypal6 <- pal_futurama(palette = c("planetexpress"), alpha = 0.7)(5)
mypal7 <- pal_startrek(palette = c("uniform"), alpha = 0.7)(5)

mycoloursP<- c(mypal, mypal2, mypal3, mypal4, mypal5, mypal6, mypal7)

#Read in
dat <- read.csv(here("data", 
                     "validation_data",
                     "20221011_RNAscope_PDGFRA_EBF1_AF.csv"))

#select detected cells only
dat <- subset(dat, dat$Name =="PathCellObject")

#select only the columns of interest

t <- select(dat, 
            Image, 
            Parent, 
            dummy_id,
            tissue,
            age_group,
            Subcellular_Channel2_Num_spots_estimated, 
            Subcellular_Channel_3_Num_spots_estimated )

#rename columns
colnames(t) [colnames(t)== "Subcellular_Channel2_Num_spots_estimated"] <- "EBF1_dots"
colnames(t) [colnames(t)== "Subcellular_Channel_3_Num_spots_estimated"] <- "PDGFRA_dots"
t$EBF1_dots <- as.numeric(t$EBF1_dots)
t$PDGFRA_dots <- as.numeric(t$PDGFRA_dots)



#select columns again
t<- select(t, Parent, dummy_id, age_group, tissue, EBF1_dots, PDGFRA_dots)

#save the reduced dataset
#write.csv(t, file = "M:/PhD/Year 4/Validation/ISH/ISH_PDGFRA_EBF1_df.csv")
#t <- read.csv("M:/PhD/Year 4/Validation/ISH/ISH_PDGFRA_EBF1_df.csv")
#Assign unique ID for each FOV - check that this works properly 
split <- split(t, t$dummy_id)
#for (i in 1:length(names(split))) {
#  temp <-split[[i]]
#  temp$FOV <- paste(names(split)[i], temp$ROI)
#  if(i == 1){
#    t <- temp
#  }else{
#    t <- rbind(temp, t)}
#}

#Add area for calculating density
t$ROI_a_um2 <- "250000"




#add logical columns 
t<- mutate(t, EBF1 = t$EBF1_dots >=1)
t<- mutate(t, PDGFRA = t$PDGFRA_dots >=3)
t<- mutate(t, dPos= t$EBF1 == TRUE & t$PDGFRA == TRUE)

#sum up positive/double positive cells per FOV (== Parent)
EBF1 <- t %>% group_by(Parent) %>%
  summarise(tot_EBF1 = sum(EBF1)) 
PDGFRA <- t %>% group_by(Parent) %>% summarise(tot_PDGFRA = sum(PDGFRA)) 
dPos <-  t %>% group_by(Parent) %>%
  summarise(tot_dPos = sum(dPos)) 
tot_nl <- t %>% group_by(Parent) %>%
  summarise(tot_nl = length(Parent))

#create df with counts per FOV
df<- merge(PDGFRA, EBF1, by = "Parent", all = TRUE)        
df<- merge(df,dPos, by = "Parent", all = TRUE) 
df <- merge(df, tot_nl, by = "Parent", all = TRUE)
t2 <- merge(df, t, by = "Parent", all = TRUE)

t2<- t2[!duplicated(t2$Parent),]


#express double pos as proportion of PDGFRA-pos
t2$prop_dPos <- t2$tot_dPos/t2$tot_PDGFRA
#calculate density 
t2$dP_density <- t2$tot_dPos/as.numeric(t2$ROI_a_um2)
t2$EBF1_density <- t2$tot_EBF1/as.numeric(t2$ROI_a_um2)
t2$PDGFRA_density <- t2$tot_PDGFRA/as.numeric(t2$ROI_a_um2)

###UPDATED##
#ANOVA works for looking at tissue and condition - data is normal so yaay
shapiro.test(t2$prop_dPos)

lin_mod <- lm(prop_dPos ~ age_group + tissue, data = t2)
summary(lin_mod)
res.aov2 <- aov(prop_dPos ~ age_group+ tissue, data = t2)
summary(res.aov2)
tukey.test <- TukeyHSD(res.aov2)

t2$age_group <- factor(t2$age_group, levels = c("young", "old"))
t2$percentage <- t2$prop_dPos *100

ggplot(t2, aes(x = age_group, y = percentage)) +
  geom_boxplot()+ 
  geom_jitter(aes(color = dummy_id),
                                 position=position_dodge(width=0.75)) +
  theme_minimal() + 
  ylab("% EBF1+ PDGFRA+ of PDGFRA+") + 
  xlab("Age Group") + 
  theme(axis.text = element_text(size = 15))  + 
  theme(legend.text = element_text(size = 12))  + 
  theme(legend.title = element_text(size = 12)) +
  scale_color_manual(values = mycoloursP)

ggplot(t2, aes(x = age_group, y = percentage)) +
  geom_boxplot()+ 
  geom_jitter(aes(color = dummy_id, shape = tissue),
              position=position_dodge(width=0.75)) +
  theme_minimal() + 
  ylab("% EBF1+ PDGFRA+ of PDGFRA+") + 
  xlab("Age Group") + 
  theme(axis.text = element_text(size = 15))  + 
  theme(legend.text = element_text(size = 12))  + 
  theme(legend.title = element_text(size = 12)) +
  scale_color_manual(values = mycoloursP)

ggplot(t2, aes(x = age_group, y = percentage)) +
  geom_boxplot()+ 
  geom_jitter(aes(shape = tissue),
              position=position_dodge(width=0.75)) +
  theme_minimal() + 
  ylab("% EBF1+ PDGFRA+ of PDGFRA+") + 
  xlab("Age Group") + 
  theme(axis.text = element_text(size = 15))  + 
  theme(legend.text = element_text(size = 12))  + 
  theme(legend.title = element_text(size = 12)) +
  scale_color_manual(values = mycoloursP)


ggplot(t2, aes(x = age_group, y = percentage)) +
  geom_boxplot()+ 
  geom_jitter(aes(color = tissue),
              position=position_dodge(width=0.75)) +
  theme_minimal() + 
  ylab("% EBF1+ PDGFRA+ of PDGFRA+") + 
  xlab("Age Group") + 
  theme(axis.text = element_text(size = 15))  + 
  theme(legend.text = element_text(size = 12))  + 
  theme(legend.title = element_text(size = 12)) +
  scale_color_manual(values = c(mycoloursP[8],mycoloursP[9]) )



# Repeat, but remove auto fluorescence signal this time
# Below I am harsh and remove any cell that has an autofluorescent signal 
# regardless of co-localisation with EBF1 signal or not.



#select only the columns of interest

t <- select(dat, 
            Image, 
            Parent, 
            dummy_id,
            age_group,
            tissue,
            Subcellular_Channel1_Num_spots_estimated,
            Subcellular_Channel2_Num_spots_estimated, 
            Subcellular_Channel_3_Num_spots_estimated )

#rename columns
colnames(t) [colnames(t)== "Subcellular_Channel1_Num_spots_estimated"] <- "AF_dots"
colnames(t) [colnames(t)== "Subcellular_Channel2_Num_spots_estimated"] <- "EBF1_dots"
colnames(t) [colnames(t)== "Subcellular_Channel_3_Num_spots_estimated"] <- "PDGFRA_dots"
t$AF_dots <- as.numeric(t$AF_dots)
t$EBF1_dots <- as.numeric(t$EBF1_dots)
t$PDGFRA_dots <- as.numeric(t$PDGFRA_dots)


#select columns again
t<- select(t, Parent, dummy_id, age_group, tissue, EBF1_dots, PDGFRA_dots, AF_dots)

#save the reduced dataset
#write.csv(t, file = "M:/PhD/Year 4/Validation/ISH/ISH_PDGFRA_EBF1_df.csv")
#t <- read.csv("M:/PhD/Year 4/Validation/ISH/ISH_PDGFRA_EBF1_df.csv")
#Assign unique ID for each FOV - check that this works properly 
split <- split(t, t$dummy_id)


#Add area for calculating density
t$ROI_a_um2 <- "250000"




#add logical columns 
t<- mutate(t, AF = t$AF_dots >=2)
t<- mutate(t, EBF1 = t$EBF1_dots >=1)
t<- mutate(t, PDGFRA = t$PDGFRA_dots >=3)
t<- mutate(t, dPos= t$EBF1 == TRUE & t$PDGFRA == TRUE)

t<- mutate(t, tripPos= t$EBF1 == TRUE & t$PDGFRA == TRUE & t$AF == TRUE)

bool <- t$tripPos == "FALSE"

summary(bool)
subs_t <- t[bool,]

#sum up positive/double positive cells per FOV (== Parent)
EBF1 <- subs_t %>% group_by(Parent) %>%
  summarise(tot_EBF1 = sum(EBF1)) 
PDGFRA <- subs_t %>% group_by(Parent) %>% summarise(tot_PDGFRA = sum(PDGFRA)) 
dPos <-  subs_t %>% group_by(Parent) %>%
  summarise(tot_dPos = sum(dPos)) 
tot_nl <- subs_t %>% group_by(Parent) %>%
  summarise(tot_nl = length(Parent))

#create df with counts per FOV
df<- merge(PDGFRA, EBF1, by = "Parent", all = TRUE)        
df<- merge(df,dPos, by = "Parent", all = TRUE) 
df <- merge(df, tot_nl, by = "Parent", all = TRUE)
t2 <- merge(df, t, by = "Parent", all = TRUE)

t2<- t2[!duplicated(t2$Parent),]


#express double pos as proportion of PDGFRA-pos
t2$prop_dPos <- t2$tot_dPos/t2$tot_PDGFRA
#calculate density 
t2$dP_density <- t2$tot_dPos/as.numeric(t2$ROI_a_um2)
t2$EBF1_density <- t2$tot_EBF1/as.numeric(t2$ROI_a_um2)
t2$PDGFRA_density <- t2$tot_PDGFRA/as.numeric(t2$ROI_a_um2)

###UPDATED##
#ANOVA works for looking at tissue and condition - data is normal so yaay
shapiro.test(t2$prop_dPos)

lin_mod <- lm(prop_dPos ~ age_group + tissue, data = t2)
summary(lin_mod)
res.aov2 <- aov(prop_dPos ~ age_group + tissue, data = t2)
summary(res.aov2)
tukey.test <- TukeyHSD(res.aov2)


lin_mod <- lm(prop_dPos ~ age_group * tissue, data = t2)
summary(lin_mod)
res.aov2 <- aov(prop_dPos ~ age_group * tissue, data = t2)
summary(res.aov2)
tukey.test <- TukeyHSD(res.aov2)

t2$age_group <- factor(t2$age_group, levels = c("young", "old"))
t2$percentage <- t2$prop_dPos *100

ggplot(t2, aes(x = age_group, y = percentage)) +
  geom_boxplot()+ 
  geom_jitter(aes(colour = dummy_id),
              position=position_dodge(width=0.75)) +
  theme_minimal() + 
  ylab("% EBF1+ PDGFRA+ of PDGFRA+") + 
  xlab("Age Group") + 
  theme(axis.text = element_text(size = 15))  + 
  theme(legend.text = element_text(size = 12))  + 
  theme(legend.title = element_text(size = 12))+
  scale_color_manual(values = mycoloursP)

ggplot(t2, aes(x = age_group, y = percentage)) +
  geom_boxplot()+ 
  geom_jitter(aes(colour = dummy_id, shape = tissue),
              position=position_dodge(width=0.75)) +
  theme_minimal() + 
  ylab("% EBF1+ PDGFRA+ of PDGFRA+") + 
  xlab("Age Group") + 
  theme(axis.text = element_text(size = 15))  + 
  theme(legend.text = element_text(size = 12))  + 
  theme(legend.title = element_text(size = 12))+
  scale_color_manual(values = mycoloursP)

ggplot(t2, aes(x = age_group, y = percentage)) +
  geom_boxplot()+ 
  geom_jitter(aes( shape = tissue),
              position=position_dodge(width=0.75)) +
  theme_minimal() + 
  ylab("% EBF1+ PDGFRA+ of PDGFRA+") + 
  xlab("Age Group") + 
  theme(axis.text = element_text(size = 15))  + 
  theme(legend.text = element_text(size = 12))  + 
  theme(legend.title = element_text(size = 12))+
  scale_color_manual(values = mycoloursP)

# Removal of cells that contain auto-fluorescent signal makes the data
# distribution non-normal and also reduces the effect size 

wilcox.test(prop_dPos ~ age_group, data = t2)



#####
####
#### CB only
cb_t2 <- subset(t2, t2$tissue == "CB")

cb_t2$prop_dPos <- cb_t2$tot_dPos/cb_t2$tot_PDGFRA
#calculate density 
cb_t2$dP_density <- cb_t2$tot_dPos/as.numeric(cb_t2$ROI_a_um2)
cb_t2$EBF1_density <- cb_t2$tot_EBF1/as.numeric(cb_t2$ROI_a_um2)
cb_t2$PDGFRA_density <- cb_t2$tot_PDGFRA/as.numeric(cb_t2$ROI_a_um2)

###UPDATED##
#ANOVA works for looking at tissue and condition - data is normal so yaay
shapiro.test(cb_t2$prop_dPos)

lin_mod <- lm(prop_dPos ~ age_group, data = cb_t2)
summary(lin_mod)
res.aov2 <- aov(prop_dPos ~ age_group, data = t2)
summary(res.aov2)
tukey.test <- TukeyHSD(res.aov2)


cb_t2$age_group <- factor(cb_t2$age_group, levels = c("young", "old"))
cb_t2$percentage <- cb_t2$prop_dPos *100

ggplot(cb_t2, aes(x = age_group, y = percentage)) +
  geom_boxplot()+ 
  geom_jitter(aes(colour = dummy_id),
              position=position_dodge(width=0.75)) +
  theme_minimal() + 
  ylab("% EBF1+ PDGFRA+ of PDGFRA+") + 
  xlab("Age Group") + 
  theme(axis.text = element_text(size = 15))  + 
  theme(legend.text = element_text(size = 12))  + 
  theme(legend.title = element_text(size = 12))+
  scale_color_manual(values = mycoloursP)


#### CSC only
csc_t2 <- subset(t2, t2$tissue == "CSC")

csc_t2$prop_dPos <- csc_t2$tot_dPos/csc_t2$tot_PDGFRA
#calculate density 
csc_t2$dP_density <- csc_t2$tot_dPos/as.numeric(csc_t2$ROI_a_um2)
csc_t2$EBF1_density <- csc_t2$tot_EBF1/as.numeric(csc_t2$ROI_a_um2)
csc_t2$PDGFRA_density <- csc_t2$tot_PDGFRA/as.numeric(csc_t2$ROI_a_um2)


shapiro.test(csc_t2$prop_dPos)

lin_mod <- lm(prop_dPos ~ age_group, data = csc_t2)
summary(lin_mod)
res.aov2 <- aov(prop_dPos ~ age_group, data = t2)
summary(res.aov2)
tukey.test <- TukeyHSD(res.aov2)


csc_t2$age_group <- factor(csc_t2$age_group, levels = c("young", "old"))
csc_t2$percentage <- csc_t2$prop_dPos *100

ggplot(csc_t2, aes(x = age_group, y = percentage)) +
  geom_boxplot()+ 
  geom_jitter(aes(colour = dummy_id),
              position=position_dodge(width=0.75)) +
  theme_minimal() + 
  ylab("% EBF1+ PDGFRA+ of PDGFRA+") + 
  xlab("Age Group") + 
  theme(axis.text = element_text(size = 15))  + 
  theme(legend.text = element_text(size = 12))  + 
  theme(legend.title = element_text(size = 12))+
  scale_color_manual(values = mycoloursP)




###### Test other autofluorescence markers to correct for 




t <- select(dat, 
            Image, 
            Parent, 
            dummy_id,
            tissue,
            age_group,
            Subcellular_Channel1_Num_spots_estimated,
            Subcellular_Channel2_Num_spots_estimated, 
            Subcellular_Channel_3_Num_spots_estimated,
            Nucleus_Alexa_488_narrow_sum,
            Nucleus_Alexa_488_narrow_std_dev,
            Nucleus_Alexa_488_narrow_max,
            Nucleus_Alexa_488_narrow_min,
            Nucleus_Alexa_488_narrow_range)




#rename columns
colnames(t) [colnames(t)== "Subcellular_Channel1_Num_spots_estimated"] <- "AF_dots"
colnames(t) [colnames(t)== "Subcellular_Channel2_Num_spots_estimated"] <- "EBF1_dots"
colnames(t) [colnames(t)== "Subcellular_Channel_3_Num_spots_estimated"] <- "PDGFRA_dots"
colnames(t)[colnames(t)== "Nucleus_Alexa_488_narrow_sum"] <- "nuc_AF_sum"
colnames(t)[colnames(t)== "Nucleus_Alexa_488_narrow_std_dev"] <- "nuc_AF_sd"
colnames(t)[colnames(t)== "Nucleus_Alexa_488_narrow_max"] <- "nuc_AF_max"
colnames(t)[colnames(t)== "Nucleus_Alexa_488_narrow_min"] <- "nuc_AF_min"
colnames(t)[colnames(t)== "Nucleus_Alexa_488_narrow_range"] <- "nuc_AF_range"

t$EBF1_dots <- as.numeric(t$EBF1_dots)
t$PDGFRA_dots <- as.numeric(t$PDGFRA_dots)
t$AF_dots <- as.numeric(t$AF_dots)
t$nuc_AF_sum <- as.numeric(t$nuc_AF_sum)
t$nuc_AF_sd <- as.numeric(t$nuc_AF_sd)
t$nuc_AF_max <- as.numeric(t$nuc_AF_max)
t$nuc_AF_min <- as.numeric(t$nuc_AF_min)
t$nuc_AF_range <- as.numeric(t$nuc_AF_range)



ggplot(t, aes(x = age_group, y = EBF1_dots))+
  geom_boxplot()+
  geom_jitter()

ggplot(t, aes(x = age_group, y = PDGFRA_dots))+
  geom_boxplot()+
  geom_jitter()

ggplot(t, aes(x = age_group, y = AF_dots))+
  geom_boxplot()+
  geom_jitter()

ggplot(t, aes(x = age_group, y = nuc_AF_sum))+
  geom_boxplot()+
  geom_jitter()

ggplot(t, aes(x = age_group, y = nuc_AF_sd))+
  geom_boxplot()+
  geom_jitter()

ggplot(t, aes(x = age_group, y = nuc_AF_max))+
  geom_boxplot()+
  geom_jitter()

ggplot(t, aes(x = age_group, y = nuc_AF_min))+
  geom_boxplot()+
  geom_jitter()

ggplot(t, aes(x = age_group, y = nuc_AF_range))+
  geom_boxplot()+
  geom_jitter()

t$mean_AF <- t$nuc_AF_range/ t$nuc_AF_sum

ggplot(t, aes(x = age_group, y = mean_AF, fill = age_group))+
  geom_boxplot()+
  geom_jitter(size = 0.1)+
  theme_minimal()+
  ylab("Mean autofluorescence")+
  xlab("Age Group")

lmod <- lm(mean_AF ~ age_group, data = t)
summary(lmod)
#select columns again
t<- select(t, Parent, dummy_id, age_group, tissue, AF_dots, 
           EBF1_dots, PDGFRA_dots, mean_AF)

#save the reduced dataset
#write.csv(t, file = "M:/PhD/Year 4/Validation/ISH/ISH_PDGFRA_EBF1_df.csv")
#t <- read.csv("M:/PhD/Year 4/Validation/ISH/ISH_PDGFRA_EBF1_df.csv")
#Assign unique ID for each FOV - check that this works properly 
split <- split(t, t$dummy_id)

#Add area for calculating density
t$ROI_a_um2 <- "250000"




#add logical columns 



t<- mutate(t, AF = t$AF_dots >=1)
t<- mutate(t, EBF1 = t$EBF1_dots >=1)
t<- mutate(t, PDGFRA = t$PDGFRA_dots >=3)
t<- mutate(t, dPos= t$EBF1 == TRUE & t$PDGFRA == TRUE)

t<- mutate(t, tripPos= t$EBF1 == TRUE & t$PDGFRA == TRUE & t$AF == TRUE)

#bool <- t$tripPos == "FALSE"

#summary(bool)
#subs_t <- t[bool,]

#sum up positive/double positive cells per FOV (== Parent)
EBF1 <- t %>% group_by(Parent) %>%
  summarise(tot_EBF1 = sum(EBF1)) 
PDGFRA <- t %>% group_by(Parent) %>% summarise(tot_PDGFRA = sum(PDGFRA)) 
dPos <-  t %>% group_by(Parent) %>%
  summarise(tot_dPos = sum(dPos)) 
tot_nl <- t %>% group_by(Parent) %>%
  summarise(tot_nl = length(Parent))

mean_af <- t %>% group_by(Parent) %>%
  summarise(tot_mean_af = mean(mean_AF))

#create df with counts per FOV
df<- merge(PDGFRA, EBF1, by = "Parent", all = TRUE)        
df<- merge(df,dPos, by = "Parent", all = TRUE) 
df <- merge(df, tot_nl, by = "Parent", all = TRUE)
df <- merge(df, t, by = "Parent", all = TRUE)
t2 <- merge(df, mean_af, by = "Parent")

t2<- t2[!duplicated(t2$Parent),]


#express double pos as proportion of PDGFRA-pos
t2$prop_dPos <- t2$tot_dPos/t2$tot_PDGFRA
#calculate density 
t2$dP_density <- t2$tot_dPos/as.numeric(t2$ROI_a_um2)
t2$EBF1_density <- t2$tot_EBF1/as.numeric(t2$ROI_a_um2)
t2$PDGFRA_density <- t2$tot_PDGFRA/as.numeric(t2$ROI_a_um2)

###UPDATED##
#ANOVA works for looking at tissue and condition - data is normal so yaay
shapiro.test(t2$prop_dPos)

lin_mod <- lm(prop_dPos ~ age_group + tissue + tot_mean_af, data = t2)
summary(lin_mod)
res.aov2 <- aov(prop_dPos ~ age_group + tissue + tot_mean_af, data = t2)
summary(res.aov2)
tukey.test <- TukeyHSD(res.aov2)




lin_mod <- lm(prop_dPos ~ age_group * tissue+ tot_mean_af, data = t2)
summary(lin_mod)
res.aov2 <- aov(prop_dPos ~ age_group * tissue+ tot_mean_af, data = t2)
summary(res.aov2)
tukey.test <- TukeyHSD(res.aov2)

t2$age_group <- factor(t2$age_group, levels = c("young", "old"))
t2$percentage <- t2$prop_dPos *100

ggplot(t2, aes(x = age_group, y = percentage)) +
  geom_boxplot()+ 
  geom_jitter(aes(colour = dummy_id),
              position=position_dodge(width=0.75)) +
  theme_minimal() + 
  ylab("% EBF1+ PDGFRA+ of PDGFRA+") + 
  xlab("Age Group") + 
  theme(axis.text = element_text(size = 15))  + 
  theme(legend.text = element_text(size = 12))  + 
  theme(legend.title = element_text(size = 12))+
  scale_color_manual(values = mycoloursP)

ggplot(t2, aes(x = age_group, y = percentage)) +
  geom_boxplot()+ 
  geom_jitter(aes(colour = dummy_id, shape = tissue),
              position=position_dodge(width=0.75)) +
  theme_minimal() + 
  ylab("% EBF1+ PDGFRA+ of PDGFRA+") + 
  xlab("Age Group") + 
  theme(axis.text = element_text(size = 15))  + 
  theme(legend.text = element_text(size = 12))  + 
  theme(legend.title = element_text(size = 12))+
  scale_color_manual(values = mycoloursP)

ggplot(t2, aes(x = age_group, y = percentage)) +
  geom_boxplot()+ 
  geom_jitter(aes( shape = tissue),
              position=position_dodge(width=0.75)) +
  theme_minimal() + 
  ylab("% EBF1+ PDGFRA+ of PDGFRA+") + 
  xlab("Age Group") + 
  theme(axis.text = element_text(size = 15))  + 
  theme(legend.text = element_text(size = 12))  + 
  theme(legend.title = element_text(size = 12))+
  scale_color_manual(values = mycoloursP)

