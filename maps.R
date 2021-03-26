# -------------------------------------------------------------------
# Load libraries
# -------------------------------------------------------------------

library(raster)
library(tmap)
library(colorRamps)
library(rasterVis)
library(grid)

library(tidyverse)
library(rgdal)
library(proj4)
library(rgeos)
library(FedData)

library(maps)
library(GISTools)
library(RColorBrewer)

# -------------------------------------------------------------------
# set directory
# -------------------------------------------------------------------

data<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting,northing,dominant.surface.rock.type)

# -------------------------------------------------------------------
# Transform spatial data (easting, northing) -> (lat, lon)
# -------------------------------------------------------------------
proj4string <- "+proj=utm +zone=11 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs "

# select easting and northing
locations <- data.frame(data$easting,data$northing)

# transform data
pj <- project(locations, proj4string, inverse=TRUE)
latlon <- data.frame(lat=pj$y, long=pj$x)

# bind to data frame
data <-bind_cols(data,latlon)

# -------------------------------------------------------------------
# Set spatial extent
# -------------------------------------------------------------------

extentClarkia <- readWKT("POLYGON((-118.25 35.45, -118.25 35.6, -118.25 35.9, -118.8 35.9, -118.8 35.45, -118.25 35.45))")
proj4string(extentClarkia) <- "+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"

# -------------------------------------------------------------------
# Get national elevation raster
# -------------------------------------------------------------------
setwd("~/Dropbox/projects/clarkiaScripts/data/ghcn")

# download National Elevation Database elevation data in extent
# default resolution is 1 arcsecond (res="1);
# to get 1/3 arcsecond (res="13)
ned_kern<-get_ned(template=extentClarkia, label="ned_kern", res="1", force.redo = F)
# try high resolution
# ned_kern<-get_ned(template=extentClarkia, label="ned_kern", res="13", force.redo = F)

myCol = colorRampPalette(brewer.pal(9, "Greys"))(20)

# -------------------------------------------------------------------
# Elevation colors
# -------------------------------------------------------------------
# Elevation colors were borrowed from the
# Kansas Geological Survey's Elevation Map of Kansas:
# http://www.kgs.ku.edu/General/elevatMap.html
myCol<-c("#06407F", "#317A9D", "#4ABEBB", "#40AE89", "#467B5D",
        "#3C6D4D", "#1A572E", "#034C00", "#045D03", "#6C975F", "#6B823A",
        "#88A237", "#C5D16B", "#DDE580", "#FFF6AE", "#FBCB81", "#F0B16A",
        "#F2B16D", "#D18338", "#B16F33", "#825337", "#66422A", "#4F2C0C")


# -------------------------------------------------------------------
# Get shapefiles of lakes and rivers
# -------------------------------------------------------------------
setwd("~/Dropbox/projects/clarkiaScripts/data/rivers")
shapeData <- rgdal::readOGR(dsn="CA_Lakes.shp",layer="CA_Lakes")
setwd("~/Dropbox/projects/clarkiaScripts/data/rivers/water_bodies")
rivers <- rgdal::readOGR(dsn="Water_Course.shp",layer="Water_Course")
bodies <- rgdal::readOGR(dsn="water_bodies.shp",layer="water_bodies")

setwd("~/Dropbox/projects/clarkiaScripts/data/rivers/rivers/Data")
rivers2 <- rgdal::readOGR(dsn="cdfg_100k_2003_6.shp",layer="cdfg_100k_2003_6")
kern2 <- subset(rivers2, NAME %in% c('Kern River'))

isabella <- subset(bodies, NAME %in% c('Lake Isabella'))

# I read in the shapefile, but I'm not sure how to work with that. But
# I do understand data frames, so that's what I'm converting it to.
#shapeData@data$id <- rownames(bodies@data)

isabella <- subset(bodies, NAME %in% c('Lake Isabella'))
isabella<-spTransform(isabella, proj4string(ned_kern))
rivers<-spTransform(rivers, proj4string(ned_kern))
bodies<-spTransform(bodies, proj4string(ned_kern))
kern2<-spTransform(kern2, proj4string(ned_kern))
#plot(kern2)
#plot(raster::crop(kern2,isabella),col="red",add=T)
#plot(isabella,add=TRUE)
kern_to_plot<-raster::crop(kern2,isabella)
#plot(kern_to_plot,col="red",add=T)
#kernFinal = kern2-kern_to_plot
kernFinal=gDifference(kern2, kern_to_plot)
#plot(kernFinal,col="blue",add=T)



xy <- data.frame(data$long,data$lat)

spdf <- SpatialPoints(coords = xy, 
                      proj4string = CRS("+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"))


# -------------------------------------------------------------------
# Plot map
# -------------------------------------------------------------------
# https://github.com/mtennekes/tmap/issues/164
# https://rstudio-pubs-static.s3.amazonaws.com/289387_6d2bbaf850764a38bfea1618e78a68ef.html
tm_shape(ned_kern, unit.size = 1)+      
  tm_grid(col = "black", n.x = 5, n.y = 5, lines = FALSE,
          labels.rot = c(0, 90),labels.size=.8) +
  tm_raster(palette = myCol,title="",
            style='cont',auto.palette.mapping=F) +
  # Add lake outline
  tm_shape(isabella) +
  tm_polygons(isabella='lightgray',col='lightgray',border.col='lightgray',alpha=.5) +
  # Add river outline
  tm_shape(kernFinal) +
  tm_lines(kernFinal='lightgray',col='lightgray',alpha=.5) +
  # Add populations
  tm_shape(spdf)+
  tm_dots(size=0.25,col='lightgray',border.col='black',border.lwd=1,shape=21)+
  # Add legend, compass, and scale bar
  tm_legend(legend.outside = TRUE,legend.outside.size=.1)+
  tm_compass(position=c(0.9,0.9))+
  tm_scale_bar(position=c('left','top'),width = 0.2,
               text.size=1,text.color='white',
               breaks = c(0,4,8,12)) 


# 
# par(mfrow=c(1,1), cex=1, mar=c(2.5,2.5,2.5,2.5))
# sp::plot(ned_kern)
# color=rgb(0,0,0,alpha=0.3)
# with(data, symbols(x=long, y=lat, circles=rep(1,20), inches=1/6, 
#                       pch=1, bg = color , add=TRUE))
# 
# legend.col <- function(col, lev){
#   
#   opar <- par
#   
#   n <- length(col)
#   
#   bx <- par("usr")
#   
#   box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
#               bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
#   box.cy <- c(bx[3], bx[3])
#   box.sy <- (bx[4] - bx[3]) / n
#   
#   xx <- rep(box.cx, each = 2)
#   
#   par(xpd = TRUE)
#   for(i in 1:n){
#     
#     yy <- c(box.cy[1] + (box.sy * (i - 1)),
#             box.cy[1] + (box.sy * (i)),
#             box.cy[1] + (box.sy * (i)),
#             box.cy[1] + (box.sy * (i - 1)))
#     polygon(xx, yy, col = col[i], border = col[i])
#     
#   }
#   par(new = TRUE)
#   plot(0, 0, type = "n",
#        ylim = c(min(lev), max(lev)),
#        yaxt = "n", ylab = "",
#        xaxt = "n", xlab = "",
#        frame.plot = FALSE)
#   axis(side = 4, las = 2, hadj=1.5,tick = FALSE, line = .25)
#   par <- opar
# }


sp::plot(ned_kern,col=myCol, box = FALSE, axes=FALSE)
colr.2 <-rgb(0,0,0,alpha=0.3)
with(data, symbols(x=long, y=lat, circles=rep(1,20), inches=1/16,
                      pch=1, bg = 'lightgray' , add=TRUE))


sp::plot(ned_kern, col=myCol, box = FALSE, axes=FALSE,mar=rep(0,4),
         useRaster=F)
colr.2 <-rgb(0,0,0,alpha=0.3)
with(data, symbols(x=long, y=lat, circles=rep(1,20), inches=1/16,
                      pch=1, bg = 'lightgray' , add=TRUE))

map.axes(cex.axis=0.8,mar=rep(0,4))
maps::map.scale(x=-118.78, y=35.89, ratio=FALSE, relwidth=0.2)
#north.arrow(xb=-116, yb=41, len=0.22, lab="N") 

library(grid)
rasterVis::levelplot(ned_kern, col.regions = myCol, cuts=length(myCol), margin=FALSE) +
  latticeExtra::layer(sp.points(spdf, col = "lightgray",pch=19)) +
  latticeExtra::layer({
    xs <- seq(-118.75, -118.65, by=.025)
    grid::grid.rect(x=xs, y=35.875,
              width=.025, height=.01,
              gp=gpar(fill=rep(c('white', 'black'), 2)),
              default.units='native')
    grid::grid.text(x= xs , y=35.885, seq(0, .2, by=.05),
              gp=gpar(cex=0.5, col = 'white'), rot=0, 
              default.units='native')
  })

raster::scalebar(d = 100, # distance in km
                 xy = c(extent(extentClarkia)[1]+0.01,extent(extentClarkia)[3]+0.15),
                 type = "bar", 
                 divs = 2, 
                 below = "km", 
                 lonlat = TRUE,
                 label = c(0,50,100), 
                 adj=c(0, -0.75), 
                 lwd = 2)

with(data, symbols(x=long, y=lat, circles=rep(1,20), inches=1/16,
                   pch=1, bg = 'lightgray' , add=TRUE))
+



# Inmap
par(usr=c(-216, -63, 22, 144))
rect(xleft =-126.2,ybottom = 23.8,xright = -100,ytop = 50.6,col = "white")
map('state', fill = FALSE, xlim = c(-125, -114), ylim = c(32.2, 42.5), xlab = "lon", ylab = "lat", add =T)

#map("usa", xlim=c(-126.2,-65.5), ylim=c(23.8,50.6),add=T)
#map("state", xlim=c(-126.2,-65.5), ylim=c(23.8,50.6),add=T, boundary = F, interior = T, lty=2)
map("state", region="california", fill=T, add=T)
points(-118.5, 35.7, bg = "white", pch = 21)
dev.off()




ned_kern<-get_ned(template=extentClarkia, label="ned_kern", res="13", force.redo = F)


# https://github.com/mtennekes/tmap/issues/164
# https://rstudio-pubs-static.s3.amazonaws.com/289387_6d2bbaf850764a38bfea1618e78a68ef.html
# Fishing ranges polygon
tm_shape(ned_kern, unit.size = 1, 
         xlim=c(-118.8,-118.6),ylim=c(35.45,35.5),
         raster.downsample = FALSE)+      
  tm_grid(col = "black", n.x = 2, n.y = 2, lines = FALSE,
          labels.rot = c(0, 90),labels.size=.8) +
  tm_raster(palette = myCol,title="",
            style='cont',auto.palette.mapping=F) +
  # Add land outline
  tm_shape(isabella) +
  tm_polygons(isabella='lightgray',col='lightgray',border.col='lightgray',alpha=.5) +
  # Add land outline
  tm_shape(kernFinal) +
  tm_lines(kernFinal='lightgray',col='lightgray',alpha=.5) +
  # add populations
  tm_shape(spdf)+
  tm_dots(size=0.25,col='lightgray',border.col='black',border.lwd=1,shape=21)+
  # Add legend, compass, and scale bar
  tm_legend(legend.outside = TRUE,legend.outside.size=.1)+
  tm_compass(position=c(0.9,0.9))+
  tm_scale_bar(position=c('left','top'),width = 0.15,
               text.size=1,text.color='white') 



tm_shape(ned_kern, unit.size = 1,xlim=c(-118.68,-118.6),ylim=c(35.5,35.6))+      
  tm_grid(col = "black", n.x = 2, n.y = 2, lines = FALSE,
          labels.rot = c(0, 90),labels.size=.8) +
  tm_raster(palette = myCol,title="",
            style='cont',auto.palette.mapping=F) +
  # Add land outline
  tm_shape(isabella) +
  tm_polygons(isabella='lightgray',col='lightgray',border.col='lightgray',alpha=.5) +
  # Add land outline
  tm_shape(kernFinal) +
  tm_lines(kernFinal='lightgray',col='lightgray',alpha=.5) +
  # add populations
  tm_shape(spdf)+
  tm_dots(size=0.25,col='lightgray',border.col='black',border.lwd=1,shape=21)+
  # Add legend, compass, and scale bar
  tm_legend(legend.outside = TRUE,legend.outside.size=.1)+
  tm_compass(position=c(0.9,0.9))+
  tm_scale_bar(position=c('left','top'),width = 0.15,
               text.size=1,text.color='white') 




tm_shape(ned_kern, unit.size = 1,xlim=c(-118.6,-118.48),ylim=c(35.525,35.65))+      
  tm_grid(col = "black", n.x = 2, n.y = 2, lines = FALSE,
          labels.rot = c(0, 90),labels.size=.8) +
  tm_raster(palette = myCol,title="",
            style='cont',auto.palette.mapping=F) +
  # Add land outline
  tm_shape(isabella) +
  tm_polygons(isabella='lightgray',col='lightgray',border.col='lightgray',alpha=.5) +
  # Add land outline
  tm_shape(kernFinal) +
  tm_lines(kernFinal='lightgray',col='lightgray',alpha=.5) +
  # add populations
  tm_shape(spdf)+
  tm_dots(size=0.25,col='lightgray',border.col='black',border.lwd=1,shape=21)+
  # Add legend, compass, and scale bar
  tm_legend(legend.outside = TRUE,legend.outside.size=.1)+
  tm_compass(position=c(0.9,0.9))+
  tm_scale_bar(position=c('left','top'),width = 0.15,
               text.size=1,text.color='white') 





tm_shape(ned_kern, unit.size = 1,xlim=c(-118.48,-118.4),ylim=c(35.55,35.625))+      
  tm_grid(col = "black", n.x = 2, n.y = 2, lines = FALSE,
          labels.rot = c(0, 90),labels.size=.8) +
  tm_raster(palette = myCol,title="",
            style='cont',auto.palette.mapping=F) +
  # Add land outline
  tm_shape(isabella) +
  tm_polygons(isabella='lightgray',col='lightgray',border.col='lightgray',alpha=.5) +
  # Add land outline
  tm_shape(kernFinal) +
  tm_lines(kernFinal='lightgray',col='lightgray',alpha=.5) +
  # add populations
  tm_shape(spdf)+
  tm_dots(size=0.25,col='lightgray',border.col='black',border.lwd=1,shape=21)+
  # Add legend, compass, and scale bar
  tm_legend(legend.outside = TRUE,legend.outside.size=.1)+
  tm_compass(position=c(0.9,0.9))+
  tm_scale_bar(position=c('left','top'),width = 0.15,
               text.size=1,text.color='white') 




tm_shape(ned_kern, unit.size = 1,xlim=c(-118.5,-118.4),ylim=c(35.68,35.75))+      
  tm_grid(col = "black", n.x = 2, n.y = 2, lines = FALSE,
          labels.rot = c(0, 90),labels.size=.8) +
  tm_raster(palette = myCol,title="",
            style='cont',auto.palette.mapping=F) +
  # Add land outline
  tm_shape(isabella) +
  tm_polygons(isabella='lightgray',col='lightgray',border.col='lightgray',alpha=.5) +
  # Add land outline
  tm_shape(kernFinal) +
  tm_lines(kernFinal='lightgray',col='lightgray',alpha=.5) +
  # add populations
  tm_shape(spdf)+
  tm_dots(size=0.25,col='lightgray',border.col='black',border.lwd=1,shape=21)+
  # Add legend, compass, and scale bar
  tm_legend(legend.outside = TRUE,legend.outside.size=.1)+
  tm_compass(position=c(0.9,0.9))+
  tm_scale_bar(position=c('left','top'),width = 0.15,
               text.size=1,text.color='white') 



tm_shape(ned_kern, unit.size = 1,xlim=c(-118.5,-118.4),ylim=c(35.75,35.85))+      
  tm_grid(col = "black", n.x = 2, n.y = 2, lines = FALSE,
          labels.rot = c(0, 90),labels.size=.8) +
  tm_raster(palette = myCol,title="",
            style='cont',auto.palette.mapping=F) +
  # Add land outline
  tm_shape(isabella) +
  tm_polygons(isabella='lightgray',col='lightgray',border.col='lightgray',alpha=.5) +
  # Add land outline
  tm_shape(kernFinal) +
  tm_lines(kernFinal='lightgray',col='lightgray',alpha=.5) +
  # add populations
  tm_shape(spdf)+
  tm_dots(size=0.25,col='lightgray',border.col='black',border.lwd=1,shape=21)+
  # Add legend, compass, and scale bar
  tm_legend(legend.outside = TRUE,legend.outside.size=.1)+
  tm_compass(position=c(0.9,0.9))+
  tm_scale_bar(position=c('left','top'),width = 0.15,
               text.size=1,text.color='white') 
