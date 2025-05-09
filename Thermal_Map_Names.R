#' ---
#' Step One: Rename Worldclim files
#' Matthew Watson
#' ---

setwd("D:/Temperature Files/Tmin") #download worldclime max and min monthly temperatures from 1961 to 2018. Set WD to where downloaded

#THIS SECTION RENAMES FILES FOR PROCESSING

p <- list.files(path = ".", pattern = "\\.tif$", full.names = TRUE)

file.rename(from=p, to=sub(pattern="-", replacement="_", p))

p <- list.files(path = ".", pattern = "\\.tif$", full.names = TRUE)

file.rename(from=p, to=sub(pattern="wc2.1_2.5m_", replacement="", p))

####seperate netcdf of CRU ts data (takes CRU mean temperature files and seperates them into individual TIFF files)
library(raster)

setwd("D:/Temperature Files/Tmean")

temp <- list.files(path ="D:/Temperature Files/Tmean/main", 
                  pattern = "\\.nc$",full.names = T)

i = 1
year=1961
while (i < 7){
  mp <- brick(temp[i])
  l=1
  u=12
  c = 1
  while (c < 11){
    n <- mp[[l:u]]
    nm <- paste0("tmp_",year,"_",".tif")
    writeRaster(n, nm, bylayer=TRUE)
    l=l+12
    u=u+12
    year=year+1
    c=c+1
  }
  i=i+1
}
