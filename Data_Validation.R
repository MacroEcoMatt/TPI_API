setwd("D:/PhD/Thesis/Data Paper/Analyses")

library(vroom)
library(dplyr)
library(ape)
library(phylolm)

#load phys data and get the lower and upper tnz values. create list of species####
phys_data <- read.csv("GlobTherm.csv") %>% filter(max_metric=="UTNZ" | min_metric == "LTNZ")

phys_data[c(9,552),7] <- "Parus gambeli"
synonym_list <- read.csv("Synonym_List.csv")
synonym_list <- synonym_list[,-1]
colnames(synonym_list)[1] <- "Binomial"

phys_data_iucn <- phys_data %>% filter(Binomial %in% synonym_list$IUCN_Binomial)
phys_data_syn <- phys_data %>% filter(!Binomial %in% phys_data_iucn$Binomial) %>%
  left_join(synonym_list) %>% filter(!is.na(IUCN_Binomial))%>%
  mutate(Binomial = IUCN_Binomial) %>% select(-"IUCN_Binomial")
phys_data_final <- rbind(phys_data_iucn, phys_data_syn)

#load thermal limits derrieved from ranges and filter to min and max values####
tpi_data <- vroom("Monthly_limits.csv")

TPI <- tpi_data %>% filter(Binomial %in% phys_data_final$Binomial) %>%
  group_by(Binomial) %>% summarize(Class = unique(Class),
                                   TPI_Max = max(TMax),
                                   TPI_Min = min(TMin),
                                   TPI_Max_mean = mean(TMax),
                                   TPI_Min_mean = mean(TMin))

TNZ <- phys_data_final %>% filter(Binomial %in% TPI$Binomial)%>%
  group_by(Binomial) %>% summarise(LTNZ = min(tmin),
                                   UTNZ = max(Tmax))
Therm_Data <- left_join(TPI,TNZ)
which(Therm_Data == "Otomys irroratus")
Therm_Data[334,"LTNZ"] <- 10
Therm_Data[334,"UTNZ"] <- 24

Bird_Data <- Therm_Data%>%filter(Class=="Aves")
Mam_Data <- Therm_Data%>%filter(Class=="Mammalia")

##
therm_amp <- vroom("D:/PhD/Thesis/Data Paper/All_Thermal.csv")%>%
  filter(Class=="Archelosauria" | Class=="Amphibia" | Class=="Lepidosauria")%>%
  select(4:8)

syn <- vroom ("D:/PhD/Thesis/TPI and API/Master_Synonym_List_All.csv")%>%rename(Binomial=Synonym)
therm_amp <- left_join(therm_amp,syn)%>%mutate(Binomial = ifelse(is.na(IUCN_Binomial), Binomial, IUCN_Binomial))%>%
  select(-IUCN_Binomial)

tpi <- vroom("D:/Thesis Projects/Climate Change Indicators/Niche Limits/TPI_API_Original/FINAL/Month_Limits.csv")%>%
  select(1,6:10)%>%filter(!is.na(AMax))%>%
  group_by(Class,Binomial)%>%
  summarise(TMin = min(TMin, na.rm=T),
            TMax = max(TMax, na.rm=T),
            AMin = min(AMin, na.rm=T),
            AMax = max(TMax, na.rm=T))

therm_amp <- left_join(therm_amp,tpi)%>%filter(!is.na(TMin))

unique(therm_amp$max_metric)
ctmax <- therm_amp %>% filter(max_metric=="ctmax")%>%rename(CTMax = Tmax)

unique(therm_amp$min_metric)
ctmin <- therm_amp %>% filter(min_metric=="ctmin")%>%rename(CTMin = tmin)%>%filter(!Class=="Amphibia")
LT50min <- therm_amp %>% filter(min_metric=="LT50")%>%rename(LT50Min = tmin)

Rep_Data_max <- ctmax%>%filter(Class=="Reptilia")
Amp_Data_max <- ctmax%>%filter(Class=="Amphibia")

Rep_Data_min <- ctmin%>%filter(Class=="Reptilia")
Amp_Data_min <- LT50min%>%filter(Class=="Amphibia")

##Load Phylogenies####
#phylogenies####
bird_tree <- read.tree("D:/Thesis Projects/Phylo Cor Traits/Phylogeny/BirdTree.tre")

oldname_b <- bird_tree$tip.label
newname_b <- bird_tree$tip.label
for (n in 1:length(newname_b)){
  newname_b[n] <- sub("([A-Za-z]+_[A-Za-z]+).*", "\\1", newname_b[n])
  newname_b[n] <- sub("_"," ", newname_b[n])
}
DF_names_b <- as.data.frame(cbind(oldname_b,newname_b))
bird_tree$tip.label<-DF_names_b[[2]][match(bird_tree$tip.label, DF_names_b[[1]])]


mammal_tree <- read.tree("D:/Thesis Projects/Climate Change Indicators/Phylogeny/MammalTree.newick")
oldname_m <- mammal_tree$tip.label
newname_m <- mammal_tree$tip.label
for (n in 1:length(newname_m)){
  newname_m[n] <- sub("([A-Za-z]+_[A-Za-z]+).*", "\\1", newname_m[n])
  newname_m[n] <- sub("_"," ", newname_m[n])
}
DF_names_m <- as.data.frame(cbind(oldname_m,newname_m))
mammal_tree$tip.label<-DF_names_m[[2]][match(mammal_tree$tip.label, DF_names_m[[1]])]

reptile_tree <- read.tree("D:/Thesis Projects/Climate Change Indicators/Phylogeny/ReptileTree.tre")
oldname_r <- reptile_tree$tip.label
newname_r <- reptile_tree$tip.label
for (n in 1:length(newname_r)){
  newname_r[n] <- sub("([A-Za-z]+_[A-Za-z]+).*", "\\1", newname_r[n])
  newname_r[n] <- sub("_"," ", newname_r[n])
}
DF_names_r <- as.data.frame(cbind(oldname_r,newname_r))
reptile_tree$tip.label<-DF_names_r[[2]][match(reptile_tree$tip.label, DF_names_r[[1]])]


amphibian_tree <- read.tree("D:/Thesis Projects/Climate Change Indicators/Phylogeny/AmphibianTree.tre")

oldname_a <- amphibian_tree$tip.label
newname_a <- amphibian_tree$tip.label
for (n in 1:length(newname_a)){
  newname_a[n] <- sub("([A-Za-z]+_[A-Za-z]+).*", "\\1", newname_a[n])
  newname_a[n] <- sub("_"," ", newname_a[n])
}
DF_names_a <- as.data.frame(cbind(oldname_a,newname_a))
amphibian_tree$tip.label<-DF_names_a[[2]][match(amphibian_tree$tip.label, DF_names_a[[1]])]
##
##loade range centroids####
centroids_all<- vroom("D:/Thesis Projects/Phylo Cor Traits/traits/Centroids_bird_north2.csv")%>%
  rbind(vroom("D:/Thesis Projects/Phylo Cor Traits/traits/Centroids_bird_south2.csv"))%>%dplyr::select(-5)%>%
  rbind(vroom("D:/Thesis Projects/Phylo Cor Traits/traits/Centroids_Updated.csv"))

cent_all <- vroom("D:/Thesis Projects/Phylo Cor Traits/traits/Centroids_all_elevation.csv")%>%
  dplyr::select(Binomial, Elevation, Elevation_Zone)

centroids_all <- left_join(centroids_all, cent_all)%>%
  mutate(Climate_Zone = ifelse(Range_Center_Lat < 23.5 & Range_Center_Lat > -23.5, "Tropical",
                               ifelse(Range_Center_Lat >= 23.5 & Range_Center_Lat < 40|
                                        Range_Center_Lat <= -23.5 & Range_Center_Lat > -40, "Subtropical",
                                      ifelse(Range_Center_Lat >= 40 & Range_Center_Lat < 60|
                                               Range_Center_Lat <= -40 & Range_Center_Lat > -60, "Temperate","Polar"))))%>%
  mutate(Global_Zone = paste0(Climate_Zone," ",Elevation_Zone))
#align species between phylogenies and datasets####
bird_phy <- bird_tree$tip.label

Bird_Data_phy <- Bird_Data %>% filter(Binomial %in% bird_phy)
bird_tree_phy <-drop.tip(bird_tree, bird_tree$tip.label[-na.omit(match(Bird_Data_phy$Binomial,bird_tree$tip.label))])
Bird_Data_phy <- Bird_Data_phy%>%distinct(Binomial, .keep_all = T) %>% 
  left_join(centroids_all)%>%
  tibble::column_to_rownames(var="Binomial")%>%
  mutate(Climate_Zone = as.factor(Climate_Zone))

Bird_Data_phy$Climate_Zone <- factor(Bird_Data_phy$Climate_Zone, levels = c("Tropical", "Sub Tropical", "Temperate", "Polar"))


mam_phy <- mammal_tree$tip.label

Mam_Data_phy <- Mam_Data %>% filter(Binomial %in% mam_phy)
mammal_tree_phy <-drop.tip(mammal_tree, mammal_tree$tip.label[-na.omit(match(Mam_Data_phy$Binomial,mammal_tree$tip.label))])

Mam_Data_phy <- Mam_Data_phy%>%distinct(Binomial, .keep_all = T)%>% 
  left_join(centroids_all)%>%
  tibble::column_to_rownames(var="Binomial")



rep_phy <- reptile_tree$tip.label

Rep_Data_max_phy <- Rep_Data_max %>% filter(Binomial %in% rep_phy)
Rep_Data_min_phy <- Rep_Data_min %>% filter(Binomial %in% rep_phy)
reptile_tree_max_phy <-drop.tip(reptile_tree, reptile_tree$tip.label[-na.omit(match(Rep_Data_max_phy$Binomial,reptile_tree$tip.label))])
reptile_tree_min_phy <-drop.tip(reptile_tree, reptile_tree$tip.label[-na.omit(match(Rep_Data_min_phy$Binomial,reptile_tree$tip.label))])

Rep_Data_max_phy <- Rep_Data_max_phy%>%
  distinct(Binomial, .keep_all = T) %>% 
  left_join(centroids_all)%>%
  distinct(Binomial, .keep_all = T) %>% 
  tibble::column_to_rownames(var="Binomial")

Rep_Data_min_phy <- Rep_Data_min_phy %>%
  distinct(Binomial, .keep_all = T)%>% 
  left_join(centroids_all)%>%
  distinct(Binomial, .keep_all = T)%>% 
  tibble::column_to_rownames(var="Binomial")

amp_phy <- amphibian_tree$tip.label

Amp_Data_max_phy <- Amp_Data_max %>% filter(Binomial %in% amp_phy)
Amp_Data_min_phy <- Amp_Data_min %>% filter(Binomial %in% amp_phy)

amp_tree_max_phy <-drop.tip(amphibian_tree, amphibian_tree$tip.label[-na.omit(match(Amp_Data_max_phy$Binomial,amphibian_tree$tip.label))])
amp_tree_min_phy <-drop.tip(amphibian_tree, amphibian_tree$tip.label[-na.omit(match(Amp_Data_min_phy$Binomial,amphibian_tree$tip.label))])

Amp_Data_max_phy <- Amp_Data_max_phy %>%
  distinct(Binomial, .keep_all = T)%>% 
  left_join(centroids_all)%>%
  distinct(Binomial, .keep_all = T)%>% 
  tibble::column_to_rownames(var="Binomial")

Amp_Data_min_phy <- Amp_Data_min_phy %>%
  distinct(Binomial, .keep_all = T)%>% 
  left_join(centroids_all)%>%
  distinct(Binomial, .keep_all = T)%>% 
  tibble::column_to_rownames(var="Binomial")

##

##Phylolm####
bird_lm_max <- phylolm(TPI_Max ~ UTNZ, 
                       data= Bird_Data_phy, phy=bird_tree_phy, model="lambda")

bird_lm_min <- phylolm(TPI_Min ~ LTNZ ,
                       data= Bird_Data_phy, phy=bird_tree_phy, model="lambda")
summary(bird_lm_max)
summary(bird_lm_min)

mammal_lm_max <- phylolm(TPI_Max ~ UTNZ ,
                         data= Mam_Data_phy, phy=mammal_tree_phy, model="lambda")
mammal_lm_min <- phylolm(TPI_Min ~ LTNZ ,
                         data= Mam_Data_phy, phy=mammal_tree_phy, model="lambda")

summary(mammal_lm_max)
summary(mammal_lm_min)


rep_lm_max <- phylolm(TMax ~ CTMax , 
                      data= Rep_Data_max_phy, phy=reptile_tree_max_phy, model="lambda")
rep_lm_min <- phylolm(TMin ~ CTMin , 
                      data= Rep_Data_min_phy, phy=reptile_tree_min_phy, model="lambda")

summary(rep_lm_max)
summary(rep_lm_min)

amp_lm_max <- phylolm(TMax ~ CTMax , 
                      data= Amp_Data_max_phy, phy=amp_tree_max_phy, model="lambda",lower.bound = 0.5)

amp_lm_min <- phylolm(TMin ~ LT50Min, 
                      data= Amp_Data_min_phy, phy=amp_tree_min_phy, model="lambda",lower.bound = 0.5)

summary(amp_lm_max)
summary(amp_lm_min)

hist(scale((Amp_Data_min_phy$Range_Center_Lat)^2))

##aic comparison

bird_lm_max_BM <- phylolm(TPI_Max ~ UTNZ, 
                          data= Bird_Data_phy, phy=bird_tree_phy, model="BM")
bird_lm_min_BM <- phylolm(TPI_Min ~ LTNZ, 
                          data= Bird_Data_phy, phy=bird_tree_phy, model="BM")

AIC(bird_lm_max)
AIC(bird_lm_max_BM)

AIC(bird_lm_min)
AIC(bird_lm_min_BM)


mammal_lm_max_BM <- phylolm(TPI_Max ~ UTNZ, 
                            data= Mam_Data_phy, phy=mammal_tree_phy, model="BM")
mammal_lm_min_BM <- phylolm(TPI_Min ~ LTNZ, 
                            data= Mam_Data_phy, phy=mammal_tree_phy, model="BM")
AIC(mammal_lm_max)
AIC(mammal_lm_max_BM)

AIC(mammal_lm_min)
AIC(mammal_lm_min_BM)


rep_lm_max_BM <- phylolm(TMax ~ CTMax , 
                         data= Rep_Data_max_phy, phy=reptile_tree_max_phy, model="BM")
rep_lm_min_BM <- phylolm(TMin ~ CTMin , 
                         data= Rep_Data_min_phy, phy=reptile_tree_min_phy, model="BM")
AIC(rep_lm_max)
AIC(rep_lm_max_BM)

AIC(rep_lm_min)
AIC(rep_lm_min_BM)


amp_lm_max_BM <- phylolm(TMax ~ CTMax, 
                         data= Amp_Data_max_phy, phy=amp_tree_max_phy, model="BM")

amp_lm_min_BM <- phylolm(TMin ~ LT50Min, 
                         data= Amp_Data_min_phy, phy=amp_tree_min_phy, model="BM")
AIC(amp_lm_max)
AIC(amp_lm_max_BM)

AIC(amp_lm_min)
AIC(amp_lm_min_BM)
#Figures####

ggplot(data=T_D, aes(x=UTNZ, y=TPI_Max))+
  geom_point(aes(color=Class), alpha=0.5, size=2)+
  geom_smooth(method = "lm", aes(color=Class), linewidth = 1.2)+
  geom_abline(slope = 1,intercept = 0, linewidth = 1.2)+
  labs(x="Upper Thermal Neutral Zone (°C)",
       y= "Upper Thermal Niche (°C)")+
  theme_light()+
  theme(axis.title = element_text(face="bold", size=14),
        axis.text = element_text(colour = "black", size=12),
        legend.position = "none",
        strip.text.x = element_text(size = 14, color = "black", face = "bold"),
        strip.background = element_rect(color="black",fill="tomato",size=1))+
  scale_color_manual(values = c("#40B0A6","tan3"))+
  facet_wrap(vars(Class), scales="free")
ggsave("D:/PhD/Thesis/Data Paper/Analyses/Figures/endo_Max.jpeg")


ggplot(data=T_D, aes(x=LTNZ, y=TPI_Min))+
  geom_point(aes(color=Class), alpha=0.5, size=2)+
  geom_smooth(method = "lm", aes(color=Class), linewidth = 1.3)+
  geom_abline(slope = 1,intercept = 0, linewidth = 1)+
  labs(x="Lower Thermal Neutral Zone (°C)",
       y= "Lower Thermal Niche (°C)")+
  theme_light()+
  theme(axis.title = element_text(face="bold", size=14),
        axis.text = element_text(colour = "black", size=12),
        legend.position = "none",
        strip.text.x = element_text(size = 14, color = "black", face = "bold"),
        strip.background = element_rect(color="black",fill="lightblue",size=1))+
  scale_color_manual(values = c("#40B0A6","tan3"))+
  facet_wrap(vars(Class), scales="free")
ggsave("D:/PhD/Thesis/Data Paper/Analyses/Figures/endo_Min.jpeg")


ggplot(data=ctmax, aes(x=CTMax, y=TMax))+
  geom_point(aes(color=Class), alpha=0.5, size=2)+
  geom_smooth(method = "lm", aes(color=Class), linewidth = 1.2)+
  geom_abline(slope = 1,intercept = 0, linewidth = 1.2)+
  labs(x="Critical Thermal Maximum (°C)",
       y= "Upper Thermal Niche (°C)")+
  theme_light()+
  theme(axis.title = element_text(face="bold", size=14),
        axis.text = element_text(colour = "black", size=12),
        legend.position = "none",
        strip.text.x = element_text(size = 14, color = "black", face = "bold"),
        strip.background = element_rect(color="black",fill="tomato",size=1))+
  scale_color_manual(values = c("purple","forestgreen"))+
  facet_wrap(vars(Class), scales="free")
ggsave("D:/PhD/Thesis/Data Paper/Analyses/Figures/exo_Max.jpeg")

a<-ggplot(data=ctmin, aes(x=CTMin, y=TMin))+
  geom_point(alpha=0.5, size=2, color="forestgreen")+
  geom_smooth(method = "lm", linewidth = 1.5, color="forestgreen")+
  geom_abline(slope = 1,intercept = 0)+
  labs(x="Critical Thermal Minimum (°C)",
       y= "Lower Thermal Position Index (°C)")+
  theme_light()+
  theme(axis.title = element_text(face="bold", size=14),
        axis.text = element_text(colour = "black", size=12),
        strip.text.x = element_text(size = 14, color = "black", face = "bold"),
        strip.background = element_rect(color="black",fill="lightblue",size=1))+
  facet_wrap(vars(Class), scales="free")
plot(a)

ggsave("D:/PhD/Thesis/Data Paper/Analyses/Figures/Reptile_Min.jpeg")

b<-ggplot(data=LT50min, aes(x=LT50Min, y=TMin))+
  geom_point(alpha=0.5, size=2, color="purple")+
  geom_smooth(method = "lm", linewidth = 1.5, color="purple")+
  geom_abline(slope = 1,intercept = 0)+
  labs(x="Minimum Lethal Temperature 50 (°C)",
       y= "Lower Thermal Position Index (°C)")+
  theme_light()+
  theme(axis.title = element_text(face="bold", size=14),
        axis.text = element_text(colour = "black", size=12),
        strip.text.x = element_text(size = 14, color = "black", face = "bold"),
        strip.background = element_rect(color="black",fill="lightblue",size=1))+
  facet_wrap(vars(Class), scales="free")
plot(b)
ggsave("D:/PhD/Thesis/Data Paper/Analyses/Figures/Amphibian_Min.jpeg")

library(ggpubr)
require(grid)
d <- ggarrange(b+ rremove("ylab") + rremove("xlab"), a + rremove("ylab") + rremove("xlab"),
               ncol=2,
               labels = NULL,
               align = "hv", 
               font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))
plot(d)

annotate_figure(d, left = textGrob("Lower Thermal Niche (°C)", rot = 90, vjust = 1, gp = gpar(cex = 1.3, fontface = 'bold')),
                bottom = textGrob("Minimum Lethal Temperature 50 (°C)", vjust = 0,gp = gpar(cex = 1.3, fontface = 'bold')))


ggsave("D:/PhD/Thesis/Data Paper/Analyses/Figures/exo_all_Min.jpeg")

##Breeding Bird Survey Analyses####
library(lme4)
library(lmerTest)
library(sjPlot)
library(dplyr)
library(ggplot2)
library(ggpubr)
###
#generate thermal limits with standard deviation####

library(rgdal)
library(terra)
library(progress)
library(dplyr)

car_chi <- vect("D:/Shapefiles/Birds/Northern/Breeding/Poecile_carolinensis.shp")

fish_crow <- vect("D:/Shapefiles/Birds/Northern/Breeding/Corvus_ossifragus.shp")

blue_jay <- vect("D:/Shapefiles/Birds/Northern/Breeding/Cyanocitta_cristata.shp")

range_lat <- data.frame(Range_Max_Lat = c(max(crds(car_chi)[, "y"]),
                                          max(crds(fish_crow)[, "y"]),
                                          max(crds(blue_jay)[, "y"])),
                        Range_Min_Lat = c(min(crds(car_chi)[, "y"]),
                                          min(crds(fish_crow)[, "y"]),
                                          min(crds(blue_jay)[, "y"])),
                        Species = c("Poecile carolinensis","Corvus ossifragus","Cyanocitta cristata"))

DFN <- data.frame(NULL)

Maxpath <- "E:/Coding Files/Climate Data/Temperature Files/Tmax/Baseline"
Minpath <- "E:/Coding Files/Climate Data/Temperature Files/Tmin/Min"

#month list Norther species
N_B_mnth <- list("_05.tif","_06.tif","_07.tif")
m <- c(5,6,7)

for(i in 1:length(N_B_mnth)){
  rmax <- rast(list.files(Maxpath, pattern = N_B_mnth[[i]], full.names = T))
  rmin <- rast(list.files(Minpath, pattern = N_B_mnth[[i]], full.names = T))
  
  maxmean_cc <- as.data.frame(extract(rmax, car_chi, fun=max, na.rm=T))
  maxmean_cc <- t(maxmean_cc[!is.infinite(rowSums(maxmean_cc)),])
  
  
  maxmean_fc <- as.data.frame(extract(rmax, fish_crow, fun=max, na.rm=T))
  maxmean_fc <- t(maxmean_fc[!is.infinite(rowSums(maxmean_fc)),])
  
  
  maxmean_bj <- as.data.frame(extract(rmax, blue_jay, fun=max, na.rm=T))
  maxmean_bj <- t(maxmean_bj[!is.infinite(rowSums(maxmean_bj)),])
  
  
  minmean_cc <- as.data.frame(extract(rmin, car_chi, fun=min, na.rm=T))
  minmean_cc <- t(minmean_cc[!is.infinite(rowSums(minmean_cc)),])
  
  
  minmean_fc <- as.data.frame(extract(rmin, fish_crow, fun=min, na.rm=T))
  minmean_fc <- t(minmean_fc[!is.infinite(rowSums(minmean_fc)),])
  
  
  minmean_bj <- as.data.frame(extract(rmin, blue_jay, fun=min, na.rm=T))
  minmean_bj <- t(minmean_bj[!is.infinite(rowSums(minmean_bj)),])
  
  
  DF <- data.frame(TMax = c(maxmean_cc[-1,1],maxmean_fc[-1,1],maxmean_bj[-1,1]),
                   TMin = c(minmean_cc[-1,1],minmean_fc[-1,1],minmean_bj[-1,1]),
                   Month = c(rep(m[i], n=45)),
                   Binomial = c("Poecile carolinensis","Poecile carolinensis","Poecile carolinensis",
                                "Poecile carolinensis","Poecile carolinensis","Poecile carolinensis",
                                "Poecile carolinensis","Poecile carolinensis","Poecile carolinensis",
                                "Poecile carolinensis","Poecile carolinensis","Poecile carolinensis",
                                "Poecile carolinensis","Poecile carolinensis","Poecile carolinensis",
                                "Corvus ossifragus", "Corvus ossifragus", "Corvus ossifragus", "Corvus ossifragus",
                                "Corvus ossifragus", "Corvus ossifragus", "Corvus ossifragus", "Corvus ossifragus",
                                "Corvus ossifragus", "Corvus ossifragus", "Corvus ossifragus", "Corvus ossifragus",
                                "Corvus ossifragus", "Corvus ossifragus", "Corvus ossifragus",
                                "Cyanocitta cristata","Cyanocitta cristata","Cyanocitta cristata","Cyanocitta cristata",
                                "Cyanocitta cristata","Cyanocitta cristata","Cyanocitta cristata","Cyanocitta cristata",
                                "Cyanocitta cristata","Cyanocitta cristata","Cyanocitta cristata","Cyanocitta cristata",
                                "Cyanocitta cristata","Cyanocitta cristata","Cyanocitta cristata"))
  
  DFN <- rbind(DFN,DF)
  gc()
}

DFN_ordered <- DFN[with(DFN, order(Binomial, Month)), ]

write.csv(DFN_ordered, "D:/PhD/Thesis/Data Paper/Revised Submission/figures/breeding_bird_Limts_all_months.csv")
#Breeding bird data####
tpi <- vroom("D:/Thesis Projects/Climate Change Indicators/Niche Limits/TPI_API_Original/FINAL/Month_Limits.csv")%>%
  dplyr::select(6:11)%>%rename(Month=MonthNumber)

lat_long <- vroom("D:/PhD/Thesis/Data Paper/Breeding Bird/routes.csv")%>%dplyr::select(1,2,3,6,7)
weather <- vroom("D:/PhD/Thesis/Data Paper/Breeding Bird/weather.csv")%>% 
  dplyr::select(1,2,3,4,5,6,7,8,11,12,13,18,19)
# Get the length of the original column
n_char <- nchar(weather$StartTime)
n_char2 <- nchar(weather$EndTime)

weather$StartTime_hour <-substr(weather$StartTime,1,nchar(weather$StartTime)-2)
weather$StartTime_Min <- substr(weather$StartTime, n_char - 1, n_char)
weather$EndTime_hour <-substr(weather$EndTime,1,nchar(weather$EndTime)-2)
weather$EndTime_Min <- substr(weather$EndTime, n_char2 - 1, n_char2)

weather <- weather %>%  mutate(removed = ifelse(TempScale == "F" & StartTemp==0, "YES","NO"))%>%
  filter(removed=="NO")%>%dplyr::select(-removed)

weather <- weather %>% mutate(StartTime_hour = as.numeric(StartTime_hour),
                              StartTime_Min = as.numeric(StartTime_Min),
                              EndTime_hour = as.numeric(EndTime_hour),
                              EndTime_Min = as.numeric(EndTime_Min),
                              total_time_min = ((EndTime_hour-StartTime_hour-1)*60) + EndTime_Min + (60-StartTime_Min))%>%
  dplyr::select(-14,-15,-16,-17)%>%
  mutate(StartTemp = ifelse(TempScale =="C" & StartTemp > 100, StartTemp-100, StartTemp),
         EndTemp = ifelse(TempScale =="C" & EndTemp > 100, EndTemp-100, EndTemp),
         StartTemp = ifelse(TempScale =="F" & StartTemp > 130, StartTemp/100, StartTemp),
         EndTemp = ifelse(TempScale =="F" & EndTemp > 130, EndTemp/100, EndTemp))%>%
  mutate(EndTemp = ifelse(TempScale=="C" & EndTemp > 45, (EndTemp - 32) * (5/9), EndTemp))%>%
  mutate(StartTemp = ifelse(TempScale == "C",StartTemp,(StartTemp - 32) * (5/9)),
         EndTemp = ifelse(TempScale == "C",EndTemp,(EndTemp - 32) * (5/9)))%>%dplyr::select(-TempScale)%>%
  mutate(StartTemp = ifelse(StartTemp > 40, (StartTemp - 32) * (5/9),StartTemp))


aou_code <- vroom("D:/PhD/Thesis/Data Paper/Breeding Bird/species_codes.csv")%>%
  mutate(Binomial = paste0(Genus," ",Species))%>%dplyr::select(AOU,Binomial)%>%
  mutate(AOU = paste0(0,AOU))

observations <- vroom("D:/PhD/Thesis/Data Paper/Breeding Bird/MigrantSummary.csv")%>%
  dplyr::select(1,2,3,4,5,6,7,14)%>%
  left_join(aou_code)%>%left_join(weather)%>%left_join(lat_long)%>%
  left_join(tpi)%>%
  mutate(Start_TPI = (StartTemp - TMin)/(TMax-TMin),
         End_TPI = (EndTemp - TMin)/(TMax-TMin),
         Site_ID = paste0(CountryNum,"_",StateNum,"_",Route))

##add in all bird observation data with climate then combine####
bird_files <- list.files("D:/PhD/Thesis/Data Paper/Breeding Bird/species count data", pattern = ".csv",full.names = T)

all_bird_data<-data.frame(NULL)

for(i in 1:length(bird_files)){
  
  cur_list <- vroom(bird_files[i])
  
  cur_list$SpeciesTotal <- rowSums(cur_list[,8:57])
  
  cur_list <- cur_list %>%dplyr::select(1:7,58)%>%
    mutate(Site_ID = paste0(CountryNum,"_",StateNum,"_",Route))
  
  all_bird_data <- rbind(all_bird_data, cur_list)
}

all_bird_data <- all_bird_data %>%distinct()%>%
  left_join(aou_code)%>%left_join(weather)%>%left_join(lat_long)%>%
  left_join(tpi)%>%
  mutate(Start_TPI = (StartTemp - TMin)/(TMax-TMin),
         End_TPI = (EndTemp - TMin)/(TMax-TMin))%>%
  filter(!is.na(Binomial))%>%filter(!is.na(TMax))

all_bird_data_temp <- all_bird_data%>%filter(!is.na(StartTemp))

obs_check <- observations %>% filter(!RouteDataID %in% all_bird_data_temp$RouteDataID)%>%
  relocate(Site_ID, .before=Binomial)

all_bird_data_temp_comb <- rbind(all_bird_data_temp,obs_check)%>%
  arrange(Site_ID,Year,Binomial)

years_check <- all_bird_data_temp_comb %>% group_by(Site_ID)%>%summarise(N_Year = n_distinct(Year))

all_bird_data_temp_comb <- left_join(all_bird_data_temp_comb, years_check)

vroom_write(all_bird_data_temp_comb, "D:/PhD/Thesis/Data Paper/Breeding Bird/All_Bird_Data_SP_Count_Final_F.csv")
###analyses####

all_bird <- vroom("D:/PhD/Thesis/Data Paper/Breeding Bird/All_Bird_Data_SP_Count_Final_F.csv")

all_bird <- all_bird %>% 
  dplyr::select(-1,-12,-15,-16)%>%
  group_by(CountryNum, StateNum, Route, RPID, Year, AOU, Site_ID, Binomial, Month)%>%
  summarise_all(mean,na.rm=T)

###species comparison####

Bb_limits <- vroom("D:/PhD/Thesis/Data Paper/Revised Submission/figures/breeding_bird_Limts_all_months.csv")[,-1]%>%
  mutate(Dataset = "Range Map")

##sp1 carolina chickadee####
"Poecile carolinensis" %in% all_bird$Binomial
##
Poecile_carolinensis <- all_bird%>%filter(Year < 2000) %>% filter(Binomial == "Poecile carolinensis")

Poecile_carolinensis <- Poecile_carolinensis %>% 
  filter(!is.na(TMax))%>%
  mutate(Start_TPI = (StartTemp - TMin) / (TMax-TMin),
         End_TPI = (EndTemp - TMin) / (TMax-TMin))%>%
  filter(!Month==8)%>%filter(!Month==4)

Poecile_carolinensis_mean <- Poecile_carolinensis%>%
  group_by(Year, Month)%>%
  summarise(Mean_Start_Temp = mean(StartTemp,na.rm=T),
            Max_Start_Temp = max(StartTemp, na.rm=T),
            Min_Start_Temp = min(StartTemp, na.rm=T),
            Mean_End_Temp = mean(EndTemp,na.rm=T),
            Max_End_Temp = max(EndTemp, na.rm=T),
            Min_End_Temp = min(EndTemp, na.rm=T),
            TMin = mean(TMin, na.rm=T),
            TMax = mean(TMax, na.rm=T),
            Lat_Max = Latitude[which.max(EndTemp)],
            Lat_Min = Latitude[which.min(StartTemp)])%>%
  ungroup()%>%filter(!is.na(TMax)&is.finite(Max_End_Temp))%>%
  ungroup()%>% mutate(Mean_Start_TPI = (Mean_Start_Temp - TMin) / (TMax-TMin),
                      Max_Start_TPI = (Max_Start_Temp - TMin) / (TMax-TMin),
                      Min_Start_TPI = (Min_Start_Temp - TMin) / (TMax-TMin),
                      Mean_End_TPI = (Mean_End_Temp - TMin) / (TMax-TMin),
                      Max_End_TPI = (Max_End_Temp - TMin) / (TMax-TMin),
                      Min_End_TPI = (Min_End_Temp - TMin) / (TMax-TMin))%>%
  filter(!Month==8)%>%filter(!Month==4)%>%
  arrange(Month,Year)
#sp2 fish crow####

"Corvus ossifragus" %in% all_bird$Binomial

Corvus_ossifragus <- all_bird%>%filter(Year < 2000) %>% filter(Binomial == "Corvus ossifragus")

Corvus_ossifragus <- Corvus_ossifragus %>% 
  filter(!is.na(TMax))%>%
  mutate(Start_TPI = (StartTemp - TMin) / (TMax-TMin),
         End_TPI = (EndTemp - TMin) / (TMax-TMin))%>%
  filter(!Month==8)%>%filter(!Month==4)%>%
  mutate(Tmax_Diff_Start =  StartTemp-TMax,
         Tmax_Diff_End =  EndTemp-TMax)

Corvus_ossifragus_mean <- Corvus_ossifragus%>%
  group_by(Year, Month)%>%
  summarise(Mean_Start_Temp = mean(StartTemp,na.rm=T),
            Max_Start_Temp = max(StartTemp, na.rm=T),
            Min_Start_Temp = min(StartTemp, na.rm=T),
            Mean_End_Temp = mean(EndTemp,na.rm=T),
            Max_End_Temp = max(EndTemp, na.rm=T),
            Min_End_Temp = min(EndTemp, na.rm=T),
            TMin = mean(TMin, na.rm=T),
            TMax = mean(TMax, na.rm=T),
            Lat_Max = Latitude[which.max(EndTemp)],
            Lat_Min = Latitude[which.min(StartTemp)])%>%
  ungroup()%>%filter(!is.na(TMax)&is.finite(Max_End_Temp))%>%
  ungroup()%>% mutate(Mean_Start_TPI = (Mean_Start_Temp - TMin) / (TMax-TMin),
                      Max_Start_TPI = (Max_Start_Temp - TMin) / (TMax-TMin),
                      Min_Start_TPI = (Min_Start_Temp - TMin) / (TMax-TMin),
                      Mean_End_TPI = (Mean_End_Temp - TMin) / (TMax-TMin),
                      Max_End_TPI = (Max_End_Temp - TMin) / (TMax-TMin),
                      Min_End_TPI = (Min_End_Temp - TMin) / (TMax-TMin))%>%
  filter(!Month==8)%>%filter(!Month==4)

#sp3 Blue Jay####
"Cyanocitta cristata" %in% all_bird$Binomial

Cyanocitta_cristata <- all_bird%>%filter(Year < 2000) %>% filter(Binomial == "Cyanocitta cristata")

Cyanocitta_cristata <- Cyanocitta_cristata %>% 
  filter(!is.na(TMax))%>%
  mutate(Start_TPI = (StartTemp - TMin) / (TMax-TMin),
         End_TPI = (EndTemp - TMin) / (TMax-TMin))%>%
  filter(!Month==8)%>%filter(!Month==4)%>%
  mutate(Tmax_Diff_Start =  StartTemp-TMax,
         Tmax_Diff_End =  EndTemp-TMax)

Cyanocitta_cristata_mean <- Cyanocitta_cristata%>%
  group_by(Year, Month)%>%
  summarise(Mean_Start_Temp = mean(StartTemp,na.rm=T),
            Max_Start_Temp = max(StartTemp, na.rm=T),
            Min_Start_Temp = min(StartTemp, na.rm=T),
            Mean_End_Temp = mean(EndTemp,na.rm=T),
            Max_End_Temp = max(EndTemp, na.rm=T),
            Min_End_Temp = min(EndTemp, na.rm=T),
            TMin = mean(TMin, na.rm=T),
            TMax = mean(TMax, na.rm=T),
            Lat_Max = Latitude[which.max(EndTemp)],
            Lat_Min = Latitude[which.min(StartTemp)])%>%
  ungroup()%>%filter(!is.na(TMax)&is.finite(Max_End_Temp))%>%
  ungroup()%>% mutate(Mean_Start_TPI = (Mean_Start_Temp - TMin) / (TMax-TMin),
                      Max_Start_TPI = (Max_Start_Temp - TMin) / (TMax-TMin),
                      Min_Start_TPI = (Min_Start_Temp - TMin) / (TMax-TMin),
                      Mean_End_TPI = (Mean_End_Temp - TMin) / (TMax-TMin),
                      Max_End_TPI = (Max_End_Temp - TMin) / (TMax-TMin),
                      Min_End_TPI = (Min_End_Temp - TMin) / (TMax-TMin))%>%
  filter(!Month==8)%>%filter(!Month==4)
##combine results for 3 species####

cc <- data.frame(TMax = Poecile_carolinensis_mean$Max_End_Temp,
                 TMin = Poecile_carolinensis_mean$Max_Start_Temp,
                 Month = Poecile_carolinensis_mean$Month,
                 Binomial = rep("Poecile carolinensis", n=nrow(Poecile_carolinensis_mean)),
                 Dataset = rep("Breeding Bird", n=nrow(Poecile_carolinensis_mean)))


fc <- data.frame(TMax = Corvus_ossifragus_mean$Max_End_Temp,
                 TMin = Corvus_ossifragus_mean$Max_Start_Temp,
                 Month = Corvus_ossifragus_mean$Month,
                 Binomial = rep("Corvus ossifragus", n=nrow(Corvus_ossifragus_mean)),
                 Dataset = rep("Breeding Bird", n=nrow(Corvus_ossifragus_mean)))

bj <- data.frame(TMax = Cyanocitta_cristata_mean$Max_End_Temp,
                 TMin = Cyanocitta_cristata_mean$Max_Start_Temp,
                 Month = Cyanocitta_cristata_mean$Month,
                 Binomial = rep("Cyanocitta cristata", n=nrow(Cyanocitta_cristata_mean)),
                 Dataset = rep("Breeding Bird", n=nrow(Cyanocitta_cristata_mean)))

Bird_compare <- rbind(Bb_limits, bj, fc, cc)

cc_stats <- rbind(cc,Bb_limits%>% filter(Binomial=="Poecile carolinensis"))
fc_stats <- rbind(fc,Bb_limits%>% filter(Binomial=="Corvus ossifragus"))
bj_stats <- rbind(bj,Bb_limits%>% filter(Binomial=="Cyanocitta cristata"))

Poc_lat_data <- Poecile_carolinensis %>% group_by(Month)%>%summarise(Max_Lat = max(Latitude),
                                                                     Min_Lat = min(Latitude))%>%
  mutate(Species= rep("Poecile carolinensis",n=3))

fishcrow_lat_data <- Corvus_ossifragus %>% group_by(Month)%>%summarise(Max_Lat = max(Latitude),
                                                                       Min_Lat = min(Latitude))%>%
  mutate(Species= rep("Corvus ossifragus",n=3))

bluejay_lat_data <- Cyanocitta_cristata %>% group_by(Month)%>%summarise(Max_Lat = max(Latitude),
                                                                        Min_Lat = min(Latitude))%>%
  mutate(Species= rep("Cyanocitta cristata",n=3))

All_Lat_data <- rbind(Poc_lat_data, fishcrow_lat_data, bluejay_lat_data)%>%left_join(range_lat)%>%
  mutate(Max_Lat_Difference = Range_Max_Lat - Max_Lat,
         Min_Lat_Difference = Range_Min_Lat - Min_Lat)%>%
  mutate(Month= ifelse(Month ==5, "May", ifelse(Month==6,"June","July")))

##analyses####

#carolina chickadee
may_cc_tmax <- t.test(cc_stats %>% filter(Month==5 & Dataset=="Range Map")%>%dplyr::select(TMax),
                      cc_stats %>% filter(Month==5 & Dataset=="Breeding Bird")%>%dplyr::select(TMax),"two.sided")
june_cc_tmax <- t.test(cc_stats %>% filter(Month==6 & Dataset=="Range Map")%>%dplyr::select(TMax),
                       cc_stats %>% filter(Month==6 & Dataset=="Breeding Bird")%>%dplyr::select(TMax),"two.sided")
july_cc_tmax <- t.test(cc_stats %>% filter(Month==7 & Dataset=="Range Map")%>%dplyr::select(TMax),
                       cc_stats %>% filter(Month==7 & Dataset=="Breeding Bird")%>%dplyr::select(TMax),"two.sided")

may_cc_tmin <- t.test(cc_stats %>% filter(Month==5 & Dataset=="Range Map")%>%dplyr::select(TMin),
                      cc_stats %>% filter(Month==5 & Dataset=="Breeding Bird")%>%dplyr::select(TMin),"two.sided")
june_cc_tmin <- t.test(cc_stats %>% filter(Month==6 & Dataset=="Range Map")%>%dplyr::select(TMin),
                       cc_stats %>% filter(Month==6 & Dataset=="Breeding Bird")%>%dplyr::select(TMin),"two.sided")
july_cc_tmin <- t.test(cc_stats %>% filter(Month==7 & Dataset=="Range Map")%>%dplyr::select(TMin),
                       cc_stats %>% filter(Month==7 & Dataset=="Breeding Bird")%>%dplyr::select(TMin),"two.sided")

#Fish Crow
may_fc_tmax <- t.test(fc_stats %>% filter(Month==5 & Dataset=="Range Map")%>%dplyr::select(TMax),
                      fc_stats %>% filter(Month==5 & Dataset=="Breeding Bird")%>%dplyr::select(TMax),"two.sided")
june_fc_tmax <- t.test(fc_stats %>% filter(Month==6 & Dataset=="Range Map")%>%dplyr::select(TMax),
                       fc_stats %>% filter(Month==6 & Dataset=="Breeding Bird")%>%dplyr::select(TMax),"two.sided")
july_fc_tmax <- t.test(fc_stats %>% filter(Month==7 & Dataset=="Range Map")%>%dplyr::select(TMax),
                       fc_stats %>% filter(Month==7 & Dataset=="Breeding Bird")%>%dplyr::select(TMax),"two.sided")

may_fc_tmin <- t.test(fc_stats %>% filter(Month==5 & Dataset=="Range Map")%>%dplyr::select(TMin),
                      fc_stats %>% filter(Month==5 & Dataset=="Breeding Bird")%>%dplyr::select(TMin),"two.sided")
june_fc_tmin <- t.test(fc_stats %>% filter(Month==6 & Dataset=="Range Map")%>%dplyr::select(TMin),
                       fc_stats %>% filter(Month==6 & Dataset=="Breeding Bird")%>%dplyr::select(TMin),"two.sided")
july_fc_tmin <- t.test(fc_stats %>% filter(Month==7 & Dataset=="Range Map")%>%dplyr::select(TMin),
                       fc_stats %>% filter(Month==7 & Dataset=="Breeding Bird")%>%dplyr::select(TMin),"two.sided")

#Blue Jay
may_bj_tmax <- t.test(bj_stats %>% filter(Month==5 & Dataset=="Range Map")%>%dplyr::select(TMax),
                      bj_stats %>% filter(Month==5 & Dataset=="Breeding Bird")%>%dplyr::select(TMax),"two.sided")
june_bj_tmax <- t.test(bj_stats %>% filter(Month==6 & Dataset=="Range Map")%>%dplyr::select(TMax),
                       bj_stats %>% filter(Month==6 & Dataset=="Breeding Bird")%>%dplyr::select(TMax),"two.sided")
july_bj_tmax <- t.test(bj_stats %>% filter(Month==7 & Dataset=="Range Map")%>%dplyr::select(TMax),
                       bj_stats %>% filter(Month==7 & Dataset=="Breeding Bird")%>%dplyr::select(TMax),"two.sided")

may_bj_tmin <- t.test(bj_stats %>% filter(Month==5 & Dataset=="Range Map")%>%dplyr::select(TMin),
                      bj_stats %>% filter(Month==5 & Dataset=="Breeding Bird")%>%dplyr::select(TMin),"two.sided")
june_bj_tmin <- t.test(bj_stats %>% filter(Month==6 & Dataset=="Range Map")%>%dplyr::select(TMin),
                       bj_stats %>% filter(Month==6 & Dataset=="Breeding Bird")%>%dplyr::select(TMin),"two.sided")
july_bj_tmin <- t.test(bj_stats %>% filter(Month==7 & Dataset=="Range Map")%>%dplyr::select(TMin),
                       bj_stats %>% filter(Month==7 & Dataset=="Breeding Bird")%>%dplyr::select(TMin),"two.sided")

#coimbine results
ttest_results <- data.frame(T_Value = c(may_cc_tmax$statistic,june_cc_tmax$statistic , july_cc_tmax$statistic, 
                                        may_cc_tmin$statistic,june_cc_tmin$statistic , july_cc_tmin$statistic,
                                        may_fc_tmax$statistic,june_fc_tmax$statistic , july_fc_tmax$statistic, 
                                        may_fc_tmin$statistic,june_fc_tmin$statistic , july_fc_tmin$statistic,
                                        may_bj_tmax$statistic,june_bj_tmax$statistic , july_bj_tmax$statistic, 
                                        may_bj_tmin$statistic,june_bj_tmin$statistic , july_bj_tmin$statistic),
                            df = c(may_cc_tmax$parameter,june_cc_tmax$parameter,july_cc_tmax$parameter,
                                   may_cc_tmin$parameter,june_cc_tmin$parameter,july_cc_tmin$parameter,
                                   may_fc_tmax$parameter,june_fc_tmax$parameter,july_fc_tmax$parameter,
                                   may_fc_tmin$parameter,june_fc_tmin$parameter,july_fc_tmin$parameter,
                                   may_bj_tmax$parameter,june_bj_tmax$parameter,july_bj_tmax$parameter,
                                   may_bj_tmin$parameter,june_bj_tmin$parameter,july_bj_tmin$parameter),
                            p_value = c(may_cc_tmax$p.value,june_cc_tmax$p.value,july_cc_tmax$p.value,
                                        may_cc_tmin$p.value,june_cc_tmin$p.value,july_cc_tmin$p.value,
                                        may_fc_tmax$p.value,june_fc_tmax$p.value,july_fc_tmax$p.value,
                                        may_fc_tmin$p.value,june_fc_tmin$p.value,july_fc_tmin$p.value,
                                        may_bj_tmax$p.value,june_bj_tmax$p.value,july_bj_tmax$p.value,
                                        may_bj_tmin$p.value,june_bj_tmin$p.value,july_bj_tmin$p.value),
                            Lower_95_CI = c(may_cc_tmax$conf.int[1],june_cc_tmax$conf.int[1],july_cc_tmax$conf.int[1],
                                            may_cc_tmin$conf.int[1],june_cc_tmin$conf.int[1],july_cc_tmin$conf.int[1],
                                            may_fc_tmax$conf.int[1],june_fc_tmax$conf.int[1],july_fc_tmax$conf.int[1],
                                            may_fc_tmin$conf.int[1],june_fc_tmin$conf.int[1],july_fc_tmin$conf.int[1],
                                            may_bj_tmax$conf.int[1],june_bj_tmax$conf.int[1],july_bj_tmax$conf.int[1],
                                            may_bj_tmin$conf.int[1],june_bj_tmin$conf.int[1],july_bj_tmin$conf.int[1]),
                            Upper_95_CI = c(may_cc_tmax$conf.int[2],june_cc_tmax$conf.int[2],july_cc_tmax$conf.int[2],
                                            may_cc_tmin$conf.int[2],june_cc_tmin$conf.int[2],july_cc_tmin$conf.int[2],
                                            may_fc_tmax$conf.int[2],june_fc_tmax$conf.int[2],july_fc_tmax$conf.int[2],
                                            may_fc_tmin$conf.int[2],june_fc_tmin$conf.int[2],july_fc_tmin$conf.int[2],
                                            may_bj_tmax$conf.int[2],june_bj_tmax$conf.int[2],july_bj_tmax$conf.int[2],
                                            may_bj_tmin$conf.int[2],june_bj_tmin$conf.int[2],july_bj_tmin$conf.int[2]),
                            Mean_Range_Map = c(may_cc_tmax$estimate[1],june_cc_tmax$estimate[1],july_cc_tmax$estimate[1],
                                               may_cc_tmin$estimate[1],june_cc_tmin$estimate[1],july_cc_tmin$estimate[1],
                                               may_fc_tmax$estimate[1],june_fc_tmax$estimate[1],july_fc_tmax$estimate[1],
                                               may_fc_tmin$estimate[1],june_fc_tmin$estimate[1],july_fc_tmin$estimate[1],
                                               may_bj_tmax$estimate[1],june_bj_tmax$estimate[1],july_bj_tmax$estimate[1],
                                               may_bj_tmin$estimate[1],june_bj_tmin$estimate[1],july_bj_tmin$estimate[1]),
                            Mean_Breeding_Bird = c(may_cc_tmax$estimate[2],june_cc_tmax$estimate[2],july_cc_tmax$estimate[2],
                                                   may_cc_tmin$estimate[2],june_cc_tmin$estimate[2],july_cc_tmin$estimate[2],
                                                   may_fc_tmax$estimate[2],june_fc_tmax$estimate[2],july_fc_tmax$estimate[2],
                                                   may_fc_tmin$estimate[2],june_fc_tmin$estimate[2],july_fc_tmin$estimate[2],
                                                   may_bj_tmax$estimate[2],june_bj_tmax$estimate[2],july_bj_tmax$estimate[2],
                                                   may_bj_tmin$estimate[2],june_bj_tmin$estimate[2],july_bj_tmin$estimate[2]),
                            Niche_Limit = c("TMax","TMax","TMax","TMin","TMin","TMin",
                                            "TMax","TMax","TMax","TMin","TMin","TMin",
                                            "TMax","TMax","TMax","TMin","TMin","TMin"),
                            Month = c("May","June","July", "May", "June", "July",
                                      "May","June","July", "May", "June", "July",
                                      "May","June","July", "May", "June", "July"),
                            Species = c("Poecile carolinensis","Poecile carolinensis","Poecile carolinensis",
                                        "Poecile carolinensis","Poecile carolinensis","Poecile carolinensis",
                                        "Corvus ossifragus","Corvus ossifragus","Corvus ossifragus",
                                        "Corvus ossifragus","Corvus ossifragus","Corvus ossifragus",
                                        "Cyanocitta cristata","Cyanocitta cristata","Cyanocitta cristata",
                                        "Cyanocitta cristata","Cyanocitta cristata","Cyanocitta cristata"))%>%
  mutate(Significant = ifelse(p_value <= 0.05, "Signficant", "Not Significant"))%>%
  left_join(All_Lat_data)
vroom_write(ttest_results, "D:/PhD/Thesis/Data Paper/Revised Submission/Range_Breeding_Bird_niche_Test.csv")

###Figures####

figure_data <- Bird_compare %>% group_by(Binomial, Month, Dataset)%>%
  summarise(TMin_sd = sd(TMin,na.rm=T),
            TMax_sd = sd(TMax, na.rm=T),
            TMin = mean(TMin,na.rm=T),
            TMax = mean(TMax, na.rm=T))%>%mutate(Month=ifelse(Month==5,"May", ifelse(Month==6, "June", "July")))%>%
  mutate(Month= as.factor(Month))


figure_data$Month <- factor(figure_data$Month, levels =  c("May", "June","July"))

##plot 1 mean and sd tmax

tmax_plot <- ggplot(figure_data, aes(x=Month, y=TMax))+
  geom_point(aes(color=Dataset), size = 2,  position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin = TMax - (2*TMax_sd), ymax = TMax + (2*TMax_sd),
                    color=Dataset),size=0.75,width=0.2, position = position_dodge(width = 0.5))+
  labs(x="Month",y="Max Temperature (Celsius)")+
  scale_color_manual(values=c("steelblue3","firebrick2","grey30"))+
  theme_minimal()+
  facet_wrap(vars(Binomial),scales="free",ncol = 3)+
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12, color="black"),
        axis.title = element_text(size = 14, color="black", face="bold"),
        legend.title =element_text(size=12,color="black",face="bold"),
        legend.text=element_text(size=12,color="black"),
        strip.text = element_text(size=14,color="black",face="italic"))
tmax_plot


##TMin plot
tmin_plot <- ggplot(figure_data, aes(x=Month, y=TMin))+
  geom_point(aes(color=Dataset), size = 2,  position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin = TMin - (2*TMin_sd), ymax = TMin + (2*TMin_sd),
                    color=Dataset),size=0.75,width=0.2, position = position_dodge(width = 0.5))+
  labs(x="Month",y="Min Temperature (Celsius)")+
  scale_color_manual(values=c("steelblue3","firebrick2","grey30"))+
  theme_minimal()+
  facet_wrap(vars(Binomial),scales="free",ncol = 3)+
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12, color="black"),
        axis.title = element_text(size = 14, color="black", face="bold"),
        legend.title =element_text(size=12,color="black",face="bold"),
        legend.text=element_text(size=12,color="black"),
        strip.text = element_text(size=14,color="black",face="italic"))
tmin_plot

##
comb_plot <- ggarrange(tmax_plot+rremove("xlab"),tmin_plot,
                       nrow=2,common.legend = T,legend = "bottom")
plot(comb_plot)
ggsave("D:/PhD/Thesis/Data Paper/Revised Submission/figures/Thermal_Niches.jpg")


##
Cyanocitta_cristata  <- Cyanocitta_cristata %>% mutate(Month_f = ifelse(Month==5,"May", ifelse(Month==6, "June", "July")))

Cyanocitta_cristata$Month_f <- factor(Cyanocitta_cristata$Month_f, levels =  c("May", "June","July"))

tpi_plot_bj <- ggplot(Cyanocitta_cristata, aes(y = End_TPI, x=Month_f))+
  geom_boxplot(aes(color=Month_f),size=1,width=0.5)+
  labs(x="Month",y="Thermal Position Index",title = "Cyanocitta cristata")+
  geom_hline(yintercept = 1, linewidth=1, linetype="dashed", color="black")+
  scale_color_manual(values=c("steelblue3","firebrick2","grey30"))+
  theme_minimal()+ylim(-0.1,1.3)+
  theme(legend.position = "none",
        axis.text = element_text(size = 14, color="black"),
        axis.title = element_text(size = 16, color="black", face="bold"),
        legend.title =element_text(size=12,color="black",face="bold"),
        legend.text=element_text(size=12,color="black"),
        plot.title = element_text(size=16,color="black", face="italic" ,hjust=0.5))

tpi_plot_bj

Poecile_carolinensis  <- Poecile_carolinensis %>% mutate(Month_f = ifelse(Month==5,"May", ifelse(Month==6, "June", "July")))
Poecile_carolinensis$Month_f <- factor(Poecile_carolinensis$Month_f, levels =  c("May", "June","July"))

tpi_plot_PC <- ggplot(Poecile_carolinensis, aes(y = End_TPI, x=Month_f))+
  
  geom_boxplot(aes(color=Month_f),size=1,width=0.5)+
  labs(x="Month",y="Thermal Position Index",title = "Poecile carolinensis")+
  geom_hline(yintercept = 1, linewidth=1, linetype="dashed", color="black")+
  scale_color_manual(values=c("steelblue3","firebrick2","grey30"))+
  theme_minimal()+ylim(-0.1,1.3)+
  theme(legend.position = "none",
        axis.text = element_text(size = 14, color="black"),
        axis.title = element_text(size = 16, color="black", face="bold"),
        legend.title =element_text(size=12,color="black",face="bold"),
        legend.text=element_text(size=12,color="black"),
        plot.title = element_text(size=16,color="black", face="italic" ,hjust=0.5))

tpi_plot_PC


Corvus_ossifragus  <- Corvus_ossifragus %>% mutate(Month_f = ifelse(Month==5,"May", ifelse(Month==6, "June", "July")))
Corvus_ossifragus$Month_f <- factor(Corvus_ossifragus$Month_f, levels =  c("May", "June","July"))

tpi_plot_fc <- ggplot(Corvus_ossifragus, aes(y = End_TPI, x=Month_f))+
  geom_boxplot(aes(color=Month_f),size=1,width=0.5)+
  labs(x="Month",y="Thermal Position Index",title = "Corvus ossifragus")+
  geom_hline(yintercept = 1, linewidth=1, linetype="dashed", color="black")+
  scale_color_manual(values=c("steelblue3","firebrick2","grey30"))+
  theme_minimal()+ylim(-0.1,1.3)+
  theme(legend.position = "none",
        axis.text = element_text(size = 14, color="black"),
        axis.title = element_text(size = 16, color="black", face="bold"),
        legend.title =element_text(size=12,color="black",face="bold"),
        legend.text=element_text(size=12,color="black"),
        plot.title = element_text(size=16,color="black", face="italic" ,hjust=0.5))

tpi_plot_fc

box_arranger <- ggarrange(tpi_plot_fc+ rremove("ylab") + rremove("xlab"),tpi_plot_bj+ rremove("ylab") + rremove("xlab"),tpi_plot_PC+ rremove("ylab") + rremove("xlab"),
                          ncol=3)

annotate_figure(box_arranger, left = textGrob("Thermal Position Index (TPI)", rot = 90, vjust = 0.5, gp = gpar(cex = 1.3, fontface = 'bold')),
                bottom = textGrob("Month", vjust = 0,gp = gpar(cex = 1.3, fontface = 'bold')))

ggsave("D:/PhD/Thesis/Data Paper/Revised Submission/figures/TPI_Boxplots.jpg")
