#' ---
#' Step 2: create aridity files
#' Matthew Watson
#' For TerraClim data this code will extract all raster layers from each
#' year and save them as geoTIFF files for future analysis
#' Ensure file paths are set to where the files should be located and set an
#' out directory to save the files
#' ---
library(terra)
library(ncdf4)
library(filenamer)

out <- "D:/Matt/Thesis Files/Universal Maps and Files/Aridity Maps/Aridity" #save location folder

#call in list of files downloaded from TerraClim (pet - potential evaporation, pre - precipitation)

pet <- list.files(path ="D:/Matt/Thesis Files/Universal Maps and Files/Aridity Maps/Evap", 
                  pattern = "\\.nc$",full.names = T)

pre <- list.files(path ="D:/Matt/Thesis Files/Universal Maps and Files/Aridity Maps/Precip/", 
                  pattern = "\\.nc$", full.names = T)

month_code <- c("Jan_","Feb_","Mar_","Apr_", "May_","June_", "July_", "Aug_","Sep_", "Oct_", "Nov_", "Dec_")
cnt=1
a = 1961
while (a < 2019) {
  evap <- nc_open(pet[cnt])
  prec <- nc_open(pre[cnt])
  
  #get variables for lat long and time
  lon <- ncvar_get(prec, "lon")
  lat <- ncvar_get(prec, "lat", verbose = F)
  
  #get variable data
  rain <- ncvar_get(prec, "ppt")
  evaporation <- ncvar_get(evap, "pet")
  
  #Fix R fill values
  fillvalue <- ncatt_get(evap, "pet", "_FillValue")
  evaporation[evaporation == fillvalue$value] <- NA
  
  fillvalue2 <- ncatt_get(prec, "ppt", "_FillValue")
  rain[rain == fillvalue2$value] <- NA
  
  b=1
  while (b < 13){
    #loop through each month in the netcdf
    p <- rain[,,b]
    e <- evaporation[,,b]
    arid <- p/e
    aridity <- raster(t(arid), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), 
                        crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
    g <- paste0("Aridity_", month_code[b], a)
    writeRaster(aridity, filename=file.path(out, as.character(a), g), "GTiff", overwrite=T)
    b = b+1
  }
  a=a+1
  cnt=cnt+1
  }
print("Done")

#averages across each year for yearly aridity maps

p <- "D:/Matt/Thesis Files/Universal Maps and Files/Aridity Maps/Aridity Yearly Average/"
y = 1961
count <- 1
while (y < 2019){
  
  evap <- rast(pet[count])
  mean_evap <- mean(evap)
  
  prec <- rast(pre[count])
  mean_pre <- mean(prec)

  arid <- mean_pre/mean_evap
  
  arid <- subst(arid, Inf, NA)
  aridyear <- paste0(p,"Aridity_",y,".tif")
  writeRaster(arid, aridyear)
  y=y+1
  count=count+1
}
