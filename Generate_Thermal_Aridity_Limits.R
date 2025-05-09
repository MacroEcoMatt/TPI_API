#' ---
#' step 4: Extract Temperature and Aridity limits
#' Matthew Watson
#' This code will extract monthly Thermal limits using worldclim data and IUCN range maps
#' Additionally a yearly average upper and lower limit will be calculated
#' due to the different range maps for birds ensure that the proper section is used
#' ---

library(terra)
library(progress)
library(dplyr)

####TPI for mammals #####
## change sp_list to the file location with the target species range maps

sp_list <- list.files("C:/Coding Files/Species Maps", pattern = "\\.shp$", full.names = T)

DF <- data.frame(Binomial= character(),
                 TMin = numeric(),
                 TMax = numeric(),
                 Month = character())

Maxpath <- "C:/Coding Files/Climate Data/MaxBase"
Minpath <- "C:/Coding Files/Climate Data/Min"

mnth <- list("_01.tif", "_02.tif", "_03.tif", "_04.tif", "_05.tif", "_06.tif",
             "_07.tif", "_08.tif", "_09.tif","_10.tif", "_11.tif", "_12.tif")

m <- list("01", "02", "03", "04", "05", "06","07", "08", "09","10", "11", "12")

###Extracts Monthly Thermal Limits

for(i in 1:length(mnth)){
  rmax <- rast(list.files(Maxpath, pattern = mnth[[i]], full.names = T))
  rmin <- rast(list.files(Minpath, pattern = mnth[[i]], full.names = T))
  pbsp <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                           total = length(sp_list),
                           complete = "=",   # Completion bar character
                           incomplete = "-", # Incomplete bar character
                           current = ">",    # Current bar character
                           clear = FALSE,    # If TRUE, clears the bar when finish
                           width = 100)
  for (sp in sp_list){
    pbsp$tick()
    polyg <- vect(sp)
    maxmean <- as.data.frame(extract(rmax, polyg, fun=max, na.rm=T))
    maxmean <- t(maxmean[is.finite(rowSums(maxmean)),])
    if (ncol(maxmean) < 2){
      mxm <- mean(maxmean[-1,])
    } else {
      mxm <- mean(colMeans(t(maxmean[-1,])))
    }
    minmean <- as.data.frame(extract(rmin, polyg, fun=min, na.rm=T))
    minmean <- t(minmean[is.finite(rowSums(minmean)),])
    if (ncol(minmean) < 2){
      mnm <- mean(minmean[-1,])
    } else {
      mnm <- mean(colMeans(t(minmean[-1,])))
    }
    allspecies <- basename(sp)
    sp_sub <- sub("_", " ", allspecies)
    bn <- sub(".shp", "", sp_sub) 
    DF[nrow(DF)+1,] = c(bn,mnm,mxm,m[[i]])
    rm(polyg, maskmax, maskmin, maxv,minv, meanmax, meanmin, allspecies, sp_sub, bn,mnm,mxm)
    gc()
  }
  rm(pbsp)
  gc()
}

## Save file ordered by species

DF_ordered <- DF[with(DF, order(Binomial, Month)), ]
write.csv(DF, "D:/Matt/Thesis Files/Universal Maps and Files/Min and Max Temp Limits/Baseline_TPI_Temperautre_Limits_Mammals.csv")

###Creates row for each species for yearly average upper and lower limit

DF_ordered$TMin <- as.numeric(DF_ordered$TMin)
DF_ordered$TMax <- as.numeric(DF_ordered$TMax)
ag <- aggregate(DF_ordered[, 2:3], list(DF_ordered$Binomial), mean)
ag<- ag %>% mutate(Month="NONE")%>% rename(Binomial = Group.1)
final <- rbind(DF_ordered,ag)
final <- final[with(final, order(Binomial, Month)),]

write.csv(final, "D:/Matt/Thesis Files/Universal Maps and Files/Min and Max Temp Limits/Baseline_TPI_Temperautre_Limits_Mammals_with_Average.csv")

#####TPI for BIRDS#####

library(rgdal)
library(terra)
library(progress)
library(dplyr)

####TPI NORTHERN#####
#Months are set to norhtern hemisphere seasons (Dec-Feb winter, Mar-May spring, June-Aug summer, Sep-Nov fall) to match range map extents

sp_list_N <- list("D:/Shapefiles/Birds/Northern/Breeding",
                 "D:/Shapefiles/Birds/Northern/Non Breeding",
                 "D:/Shapefiles/Birds/Northern/Migration Non Breeding",
                 "D:/Shapefiles/Birds/Northern/Migration Breeding")

DFN <- data.frame(Binomial= character(),
                 TMin = numeric(),
                 TMax = numeric(),
                 Month = character())

Maxpath <- "C:/Coding Files/Climate Data/MaxBase"
Minpath <- "C:/Coding Files/Climate Data/Min"

#month list Norther species
N_B_mnth <- list("_05.tif", "_06.tif","_07.tif", "_08.tif", "_09.tif")

N_NB_mnth <-list("_01.tif", "_02.tif","_12.tif")

N_NBM_mnth <-list("_03.tif","_11.tif")

N_BM_mnth <-list("_04.tif","_10.tif")

monthslistN <- (N_B_mnth, N_NB_mnth, N_NBM_mnth, N_BM_mnth)
cn=1

while (cn < 5){
  mn <- monthlistN[[cn]]
  splist <- list.files(sp_list_N[[cn]], pattern = "\\.shp$", full.name=T)
  for(i in 1:length(mn)){
    rmax <- rast(list.files(Maxpath, pattern = mn[[i]], full.names = T))
    rmin <- rast(list.files(Minpath, pattern = mn[[i]], full.names = T))
    pbsp <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                            total = length(sp_list),
                            complete = "=",   # Completion bar character
                            incomplete = "-", # Incomplete bar character
                            current = ">",    # Current bar character
                            clear = FALSE,    # If TRUE, clears the bar when finish
                            width = 100)
    for (sp in splist){
      pbsp$tick()
      polyg <- vect(sp)
      maxmean <- as.data.frame(extract(rmax, polyg, fun=max, na.rm=T))
      maxmean <- t(maxmean[!is.infinite(rowSums(maxmean)),])
      if (ncol(maxmean) < 2){
        mxm <- mean(maxmean[-1,])
      } else {
        mxm <- mean(colMeans(t(maxmean[-1,])))
      }
      minmean <- as.data.frame(extract(rmin, polyg, fun=min, na.rm=T))
      minmean <- t(minmean[!is.infinite(rowSums(minmean)),])
      if (ncol(minmean) < 2){
        mnm <- mean(minmean[-1,])
      } else {
        mnm <- mean(colMeans(t(minmean[-1,])))
      }
      allspecies <- basename(sp)
      sp_sub <- sub("_", " ", allspecies)
      bn <- sub(".shp", "", sp_sub)
      m <- sub("_","",mn[[i]])
      m <- sub(".tif","", m)
      DFN[nrow(DFN)+1,] = c(bn,mnm,mxm,m)
      rm(polyg, maskmax, maskmin, maxv,minv, meanmax, meanmin, allspecies, sp_sub, bn,mnm,mxm)
      gc()
    }
    rm(pbsp)
    gc()
  }
  cn=cn+1
}  

DFN_ordered <- DFN[with(DFN, order(Binomial, Month)), ]
write.csv(DFN_ordered, "D:/Matt/Thesis Files/Universal Maps and Files/Min and Max Temp Limits/Baseline_TPI_Temperautre_Limits_Birds_north.csv")
DFN_ordered$TMin <- as.numeric(DFN_ordered$TMin)
DFN_ordered$TMax <- as.numeric(DFN_ordered$TMax)

ag <- aggregate(DFN_ordered[, 2:3], list(DFN_ordered$Binomial), mean)
ag<- ag %>% mutate(Month="NONE")%>% rename(Binomial = Group.1)
final <- rbind(DFN_ordered,ag)
final <- final[with(final, order(Binomial, Month)),]

write.csv(final, "D:/Matt/Thesis Files/Universal Maps and Files/Min and Max Temp Limits/Baseline_TPI_Temperautre_Limits_Birds_North_with_Average.csv")

#month list Southern species
#Months are set to southern hemisphere seasons (Dec-Feb summer, Mar-May fall, June-Aug winter, Sep-Nov spring) to match range map extents

sp_list_S <- list("D:/Shapefiles/Birds/Southern/Breeding",
                 "D:/Shapefiles/Birds/Southern/Non Breeding",
                 "D:/Shapefiles/Birds/Southern/Migration Non Breeding",
                 "D:/Shapefiles/Birds/Southern/Migration Breeding")

S_B_mnth <- list("_11.tif", "_12.tif","_01.tif", "_02.tif", "_03.tif")

S_NB_mnth <-list("_06.tif", "_07.tif","_08.tif")

S_NBM_mnth <-list("_05.tif","_09.tif")

S_BM_mnth <-list("_04.tif","_10.tif")

DFS <- data.frame(Binomial= character(),
                 TMin = numeric(),
                 TMax = numeric(),
                 Month = character())

monthslistS <- (S_B_mnth, S_NB_mnth, S_NBM_mnth, S_BM_mnth)
cn=1

while (cn < 5){
  mn <- monthlistS[[cn]]
  splist <- list.files(sp_list_S[[cn]], pattern = "\\.shp$", full.name=T)
  for(i in 1:length(mn)){
    rmax <- rast(list.files(Maxpath, pattern = mn[[i]], full.names = T))
    rmin <- rast(list.files(Minpath, pattern = mn[[i]], full.names = T))
    pbsp <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                            total = length(sp_list),
                            complete = "=",   # Completion bar character
                            incomplete = "-", # Incomplete bar character
                            current = ">",    # Current bar character
                            clear = FALSE,    # If TRUE, clears the bar when finish
                            width = 100)
    for (sp in splist){
      pbsp$tick()
      polyg <- vect(sp)
      maxmean <- as.data.frame(extract(rmax, polyg, fun=max, na.rm=T))
      maxmean <- t(maxmean[!is.infinite(rowSums(maxmean)),])
      if (ncol(maxmean) < 2){
        mxm <- mean(maxmean[-1,])
      } else {
        mxm <- mean(colMeans(t(maxmean[-1,])))
      }
      minmean <- as.data.frame(extract(rmin, polyg, fun=min, na.rm=T))
      minmean <- t(minmean[!is.infinite(rowSums(minmean)),])
      if (ncol(minmean) < 2){
        mnm <- mean(minmean[-1,])
      } else {
        mnm <- mean(colMeans(t(minmean[-1,])))
      }
      allspecies <- basename(sp)
      sp_sub <- sub("_", " ", allspecies)
      bn <- sub(".shp", "", sp_sub)
      m <- sub("_","",mn[[i]])
      m <- sub(".tif","", m)
      DFS[nrow(DFS)+1,] = c(bn,mnm,mxm,m)
      rm(polyg, maskmax, maskmin, maxv,minv, meanmax, meanmin, allspecies, sp_sub, bn,mnm,mxm)
      gc()
    }
    rm(pbsp)
    gc()
  }
  cn=cn+1
}  

DFS_ordered <- DFS[with(DFS, order(Binomial, Month)), ]
write.csv(DFS_ordered, "D:/Matt/Thesis Files/Universal Maps and Files/Min and Max Temp Limits/Baseline_TPI_Temperautre_Limits_Birds_south.csv")
DFS_ordered$TMin <- as.numeric(DFS_ordered$TMin)
DFS_ordered$TMax <- as.numeric(DFS_ordered$TMax)

ags <- aggregate(DFS_ordered[, 2:3], list(DFS_ordered$Binomial), mean)
ags<- ags %>% mutate(Month="NONE")%>% rename(Binomial = Group.1)
finals <- rbind(DFS_ordered,ags)
finals <- final[with(final, order(Binomial, Month)),]

write.csv(finals, "D:/Matt/Thesis Files/Universal Maps and Files/Min and Max Temp Limits/Baseline_TPI_Temperautre_Limits_Birds_South_with_Average.csv")


library(terra)
library(progress)
library(dplyr)

####TPI for mammals and reptiles #####
## change sp_list to the file location with the target species range maps

sp_list <- list.files("E:/Coding Files/Range Maps/Mammals/Individuals", 
                      pattern = "\\.shp$", full.names = T)

DF <- data.frame(Binomial= character(),
                 AMin = numeric(),
                 AMax = numeric(),
                 Month = character())

arid_path <- "E:/Coding Files/Climate Data/Environmental Maps/Baseline/"

m <- list("Jan", "Feb", "Mar", "Apr", "May", "June","July", "Aug", "Sep","Oct", "Nov", "Dec")

for(i in 1:length(m)){
  rast_arid <- rast(list.files(arid_path, pattern = m[[i]], full.names = T))
  pbsp <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                           total = length(sp_list),
                           complete = "=",   # Completion bar character
                           incomplete = "-", # Incomplete bar character
                           current = ">",    # Current bar character
                           clear = FALSE,    # If TRUE, clears the bar when finish
                           width = 100)
  for (sp in sp_list){
    pbsp$tick()
    
    polyg <- vect(sp)
    quant <- crop(rast_arid, polyg)
    quant <- mask(rast_arid, polyg)
    
    lower_q <- sapply((global(quant, quantile, probs=c(0.25), na.rm=T)),as.numeric)
    upper_q <- sapply((global(quant, quantile, probs=c(0.75), na.rm=T)),as.numeric)
    
    if (ncol(lower_q) < 2){
      lower_q <- min(lower_q)
    } else {
      lower_q <- min(apply(lower_q, 1, min, na.rm=TRUE))
    }
    
    if (ncol(upper_q) < 2){
      upper_q <- max(upper_q)
    } else {
      upper_q <- max(apply(upper_q, 1, max, na.rm=TRUE))
    }
    
    IQRr <- as.data.frame(extract(rast_arid, polyg, fun=IQR, na.rm=T))[,-1]
    IQRr <- t(IQRr[is.finite(rowSums(IQRr)),])
    if (ncol(IQRr) < 2){
      IQRr <- max(IQRr)
    } else {
      IQRr <- max(apply(IQRr, 1, max, na.rm=TRUE))
    }
    
    
    maxair <- as.data.frame(extract(rast_arid, polyg, fun=max, na.rm=T))[,-1]
    maxair <- t(maxair[is.finite(rowSums(maxair)),])
    maxair[maxair > (upper_q + (1.5*IQRr))] <- (upper_q + (1.5*IQRr))
    if (ncol(maxair) < 2){
      maxair <- mean(maxair)
    } else {
      maxair <- mean(apply(maxair, 1, max, na.rm=TRUE))
    }
    
    
    minair <- as.data.frame(extract(rast_arid, polyg, fun=min, na.rm=T))[,-1]
    minair <- t(minair[is.finite(rowSums(minair)),])
    minair[minair < (lower_q - (1.5*IQRr))] <- (lower_q - (1.5*IQRr))
    if (ncol(minair) < 2){
      minair <- mean(minair)
    } else {
      minair <- mean(apply(minair, 1, min, na.rm=TRUE))
    }
    
    allspecies <- basename(sp)
    sp_sub <- sub("_", " ", allspecies)
    bn <- sub(".shp", "", sp_sub) 
    DF[nrow(DF)+1,] = c(bn,minair,maxair,m[[i]])
    rm(polyg, quant, minair, maxair,IQRr, lower_q, upper_q, allspecies, sp_sub, bn)
    gc()
  }
  rm(pbsp)
}

DF_ordered <- DF[with(DF, order(Binomial, Month)), ]
write.csv(DF, "C:/Users/matth/OneDrive/Documents/PhD/Thesis/TPI and API/Updated AI with 95% ci to reduce 100s/Jan_Mammals.csv")
DF_ordered$TMin <- as.numeric(DF_ordered$TMin)
DF_ordered$TMax <- as.numeric(DF_ordered$TMax)

ag <- aggregate(DF_ordered[, 2:3], list(DF_ordered$Binomial), mean)
ag<- ag %>% mutate(Month="NONE")%>% rename(Binomial = Group.1)
final <- rbind(DF_ordered,ag)
final <- final[with(final, order(Binomial, Month)),]

write.csv(final, "D:/Matt/Thesis Files/Universal Maps and Files/Min and Max Temp Limits/Baseline_TPI_Temperautre_Limits_Mammals_with_Average.csv")


###Reptiles###
sp_listr <- list.files("E:/Coding Files/Range Maps/Reptiles/All Species", 
                      pattern = "\\.shp$", full.names = T)

DFR <- data.frame(Binomial= character(),
                 AMin = numeric(),
                 AMax = numeric(),
                 Month = character())

arid_path <- "E:/Coding Files/Climate Data/Environmental Maps/Baseline/"

m <- list("Jan", "Feb", "Mar", "Apr", "May", "June","July", "Aug", "Sep","Oct", "Nov", "Dec")

for(i in 1:length(m)){
  rast_arid <- rast(list.files(arid_path, pattern = m[[i]], full.names = T))
  pbsp <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                           total = length(sp_list),
                           complete = "=",   # Completion bar character
                           incomplete = "-", # Incomplete bar character
                           current = ">",    # Current bar character
                           clear = FALSE,    # If TRUE, clears the bar when finish
                           width = 100)
  for (sp in sp_listr){
    pbsp$tick()
    
    polyg <- vect(sp)
    quant <- crop(rast_arid, polyg)
    quant <- mask(rast_arid, polyg)
    
    lower_q <- sapply((global(quant, quantile, probs=c(0.25), na.rm=T)),as.numeric)
    upper_q <- sapply((global(quant, quantile, probs=c(0.75), na.rm=T)),as.numeric)
    
    if (ncol(lower_q) < 2){
      lower_q <- min(lower_q)
    } else {
      lower_q <- min(apply(lower_q, 1, min, na.rm=TRUE))
    }
    
    if (ncol(upper_q) < 2){
      upper_q <- max(upper_q)
    } else {
      upper_q <- max(apply(upper_q, 1, max, na.rm=TRUE))
    }
    
    IQRr <- as.data.frame(extract(rast_arid, polyg, fun=IQR, na.rm=T))[,-1]
    IQRr <- t(IQRr[is.finite(rowSums(IQRr)),])
    if (ncol(IQRr) < 2){
      IQRr <- max(IQRr)
    } else {
      IQRr <- max(apply(IQRr, 1, max, na.rm=TRUE))
    }
    
    
    maxair <- as.data.frame(extract(rast_arid, polyg, fun=max, na.rm=T))[,-1]
    maxair <- t(maxair[is.finite(rowSums(maxair)),])
    maxair[maxair > (upper_q + (1.5*IQRr))] <- (upper_q + (1.5*IQRr))
    if (ncol(maxair) < 2){
      maxair <- mean(maxair)
    } else {
      maxair <- mean(apply(maxair, 1, max, na.rm=TRUE))
    }
    
    
    minair <- as.data.frame(extract(rast_arid, polyg, fun=min, na.rm=T))[,-1]
    minair <- t(minair[is.finite(rowSums(minair)),])
    minair[minair < (lower_q - (1.5*IQRr))] <- (lower_q - (1.5*IQRr))
    if (ncol(minair) < 2){
      minair <- mean(minair)
    } else {
      minair <- mean(apply(minair, 1, min, na.rm=TRUE))
    }
    
    allspecies <- basename(sp)
    sp_sub <- sub("_", " ", allspecies)
    bn <- sub(".shp", "", sp_sub) 
    DFR[nrow(DFR)+1,] = c(bn,minair,maxair,m[[i]])
    rm(polyg, quant, minair, maxair,IQRr, lower_q, upper_q, allspecies, sp_sub, bn)
    gc()
  }
  rm(pbsp)
}

DF_orderedr <- DFR[with(DFR, order(Binomial)), ]
write.csv(DF_orderedr, "C:/Users/matth/OneDrive/Documents/PhD/Thesis/TPI and API/Updated AI with 95% ci to reduce 100s/Jan_Mammals.csv")
DF_orderedr$AMin <- as.numeric(DF_orderedr$AMin)
DF_orderedr$AMax <- as.numeric(DF_orderedr$AMax)

ag <- aggregate(DF_orderedr[, 2:3], list(DF_orderedr$Binomial), mean)
ag<- ag %>% mutate(Month="NONE")%>% rename(Binomial = Group.1)
final <- rbind(DF_orderedr,ag)
final <- final[with(final, order(Binomial)),]

write.csv(final, "D:/Matt/Thesis Files/Universal Maps and Files/Min and Max Temp Limits/Baseline_TPI_Temperautre_Limits_Mammals_with_Average.csv")


###Birds North###

library(terra)
library(progress)
library(dplyr)

sp_list_N <- list("D:/Shapefiles/Birds/Northern/Breeding",
                  "D:/Shapefiles/Birds/Northern/Non Breeding",
                  "D:/Shapefiles/Birds/Northern/Non Breeding Migration",
                  "D:/Shapefiles/Birds/Northern/Breeding Migration")

DFN <- data.frame(Binomial= character(),
                  AMin = numeric(),
                  AMax = numeric(),
                  Month = character())

arid_path <- "E:/Coding Files/Climate Data/Environmental Maps/Baseline/"

#month list Norther species
N_B_mnth <- list("May", "June","July", "Aug", "Sep")

N_NB_mnth <-list("Jan", "Feb","Dec")

N_NBM_mnth <-list("Mar","Nov")

N_BM_mnth <-list("Apr","Oct")

monthslistN <- list(N_B_mnth, N_NB_mnth, N_NBM_mnth, N_BM_mnth)
cn=3

while (cn < 4){
  mn <- monthslistN[[cn]]
  splist <- list.files(sp_list_N[[cn]], pattern = "\\.shp$", full.name=T)
  for(i in 1:length(mn)){
    pbsp <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                             total = length(splist),
                             complete = "=",   # Completion bar character
                             incomplete = "-", # Incomplete bar character
                             current = ">",    # Current bar character
                             clear = FALSE,    # If TRUE, clears the bar when finish
                             width = 100)
    
    rast_arid <- rast(list.files(arid_path, pattern = mn[[i]], full.names = T))
    for (sp in splist){
      pbsp$tick()
      
      polyg <- vect(sp)
      quant <- crop(rast_arid, polyg)
      quant <- mask(quant, polyg)
      
      lower_q <- sapply((global(quant, quantile, probs=c(0.25), na.rm=T)),as.numeric)
      upper_q <- sapply((global(quant, quantile, probs=c(0.75), na.rm=T)),as.numeric)
      
      if (ncol(lower_q) < 2){
        lower_q <- min(lower_q)
      } else {
        lower_q <- min(apply(lower_q, 1, min, na.rm=TRUE))
      }
      
      if (ncol(upper_q) < 2){
        upper_q <- max(upper_q)
      } else {
        upper_q <- max(apply(upper_q, 1, max, na.rm=TRUE))
      }
      
      IQRr <- as.data.frame(extract(rast_arid, polyg, fun=IQR, na.rm=T))[,-1]
      IQRr <- t(IQRr[is.finite(rowSums(IQRr)),])
      if (ncol(IQRr) < 2){
        IQRr <- max(IQRr)
      } else {
        IQRr <- max(apply(IQRr, 1, max, na.rm=TRUE))
      }
      
      
      maxair <- as.data.frame(extract(rast_arid, polyg, fun=max, na.rm=T))[,-1]
      maxair <- t(maxair[is.finite(rowSums(maxair)),])
      maxair[maxair > (upper_q + (1.5*IQRr))] <- (upper_q + (1.5*IQRr))
      if (ncol(maxair) < 2){
        maxair <- mean(maxair)
      } else {
        maxair <- mean(apply(maxair, 1, max, na.rm=TRUE))
      }
      
      
      minair <- as.data.frame(extract(rast_arid, polyg, fun=min, na.rm=T))[,-1]
      minair <- t(minair[is.finite(rowSums(minair)),])
      if (ncol(minair) < 2){
        minair <- mean(minair)
      } else {
        minair <- mean(apply(minair, 1, min, na.rm=TRUE))
      }
      
      allspecies <- basename(sp)
      sp_sub <- sub("_", " ", allspecies)
      bn <- sub(".shp", "", sp_sub) 
      DFN[nrow(DFN)+1,] = c(bn,minair,maxair,mn[[i]])
      rm(polyg, quant, minair, maxair,IQRr, lower_q, upper_q, allspecies, sp_sub, bn)
      gc()
    }
  }
  cn=cn+1
}

DF_ordered2 <- DFN[with(DFN, order(Binomial)), ]
write.csv(DFN, "C:/Users/matth/OneDrive/Documents/PhD/Thesis/TPI and API/Updated AI with 95% ci to reduce 100s/Arid_Birds_North.csv")
DF_ordered2$TMin <- as.numeric(DF_ordered2$TMin)
DF_ordered2$TMax <- as.numeric(DF_ordered2$TMax)
colnames(DF_ordered2)<- c("Binomial","AMin","AMax","Month")

ag <- aggregate(DF_ordered2[, 2:3], list(DF_ordered$Binomial), mean)
ag<- ag %>% mutate(Month="NONE")%>% rename(Binomial = Group.1)
final <- rbind(DF_ordered2,ag)
final <- final[with(final, order(Binomial)),]

write.csv(final, "C:/Users/matth/OneDrive/Documents/PhD/Thesis/TPI and API/Updated AI with 95% ci to reduce 100s/Arid_Birds_North_Final.csv")

###Birds South###

library(terra)
library(progress)
library(dplyr)

sp_list_S <- list("D:/Shapefiles/Birds/Southern/Breeding",
                  "D:/Shapefiles/Birds/Southern/Non Breeding",
                  "D:/Shapefiles/Birds/Southern/Non Breeding Migration",
                  "D:/Shapefiles/Birds/Southern/Breeding Migration")

S_B_mnth <- list("Nov", "Dec","Jan", "Feb", "Mar")

S_NB_mnth <-list("June", "July","Aug")

S_NBM_mnth <-list("May","Sep")

S_BM_mnth <-list("Apr","Oct")

DFS <- data.frame(Binomial= character(),
                  AMax = numeric(),
                  AMax = numeric(),
                  Month = character())

arid_path <- "E:/Coding Files/Climate Data/Environmental Maps/Baseline/"

monthslistS <- list(S_B_mnth, S_NB_mnth, S_NBM_mnth, S_BM_mnth)
cn=1

while (cn < 5){
  mn <- monthslistS[[cn]]
  splist <- list.files(sp_list_S[[cn]], pattern = "\\.shp$", full.name=T)
  for(i in 1:length(mn)){
    pbsp <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                             total = length(splist),
                             complete = "=",   # Completion bar character
                             incomplete = "-", # Incomplete bar character
                             current = ">",    # Current bar character
                             clear = FALSE,    # If TRUE, clears the bar when finish
                             width = 100)
    
    rast_arid <- rast(list.files(arid_path, pattern = mn[[i]], full.names = T))
    for (sp in splist){
      pbsp$tick()
      
      polyg <- vect(sp)
      quant <- crop(rast_arid, polyg)
      quant <- mask(quant, polyg)
      
      lower_q <- sapply((global(quant, quantile, probs=c(0.25), na.rm=T)),as.numeric)
      upper_q <- sapply((global(quant, quantile, probs=c(0.75), na.rm=T)),as.numeric)
      
      if (ncol(lower_q) < 2){
        lower_q <- min(lower_q)
      } else {
        lower_q <- min(apply(lower_q, 1, min, na.rm=TRUE))
      }
      
      if (ncol(upper_q) < 2){
        upper_q <- max(upper_q)
      } else {
        upper_q <- max(apply(upper_q, 1, max, na.rm=TRUE))
      }
      
      IQRr <- as.data.frame(extract(rast_arid, polyg, fun=IQR, na.rm=T))[,-1]
      IQRr <- t(IQRr[is.finite(rowSums(IQRr)),])
      if (ncol(IQRr) < 2){
        IQRr <- max(IQRr)
      } else {
        IQRr <- max(apply(IQRr, 1, max, na.rm=TRUE))
      }
      
      
      maxair <- as.data.frame(extract(rast_arid, polyg, fun=max, na.rm=T))[,-1]
      maxair <- t(maxair[is.finite(rowSums(maxair)),])
      maxair[maxair > (upper_q + (1.5*IQRr))] <- (upper_q + (1.5*IQRr))
      if (ncol(maxair) < 2){
        maxair <- mean(maxair)
      } else {
        maxair <- mean(apply(maxair, 1, max, na.rm=TRUE))
      }
      
      
      minair <- as.data.frame(extract(rast_arid, polyg, fun=min, na.rm=T))[,-1]
      minair <- t(minair[is.finite(rowSums(minair)),])
      if (ncol(minair) < 2){
        minair <- mean(minair)
      } else {
        minair <- mean(apply(minair, 1, min, na.rm=TRUE))
      }
      
      allspecies <- basename(sp)
      sp_sub <- sub("_", " ", allspecies)
      bn <- sub(".shp", "", sp_sub) 
      DFS[nrow(DFS)+1,] = c(bn,minair,maxair,mn[[i]])
      rm(polyg, quant, minair, maxair,IQRr, lower_q, upper_q, allspecies, sp_sub, bn)
      gc()
    }
  }
  cn=cn+1
}

DFS_ordered <- DFS[with(DFS, order(Binomial)), ]
colnames(DFS_ordered)<- c("Binomial","AMin","AMax","Month")
write.csv(DFS, "C:/Users/matth/OneDrive/Documents/PhD/Thesis/TPI and API/Updated AI with 95% ci to reduce 100s/Arid_Birds_South.csv")
DFS_ordered$AMin <- as.numeric(DFS_ordered$AMin)
DFS_ordered$AMax <- as.numeric(DFS_ordered$AMax)

ag <- aggregate(DFS_ordered[, 2:3], list(DFS_ordered$Binomial), mean)
ag<- ag %>% mutate(Month="NONE")%>% rename(Binomial = Group.1)
final <- rbind(DFS_ordered,ag)
final <- final[with(final, order(Binomial)),]

write.csv(final, "C:/Users/matth/OneDrive/Documents/PhD/Thesis/TPI and API/Updated AI with 95% ci to reduce 100s/Arid_Birds_South_Final.csv")
