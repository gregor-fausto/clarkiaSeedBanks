library(GISTools)
plot.new()
#par(usr=c(-200, , 22, 144))
#rect(xleft =-200,ybottom = 23.8,xright = -100,ytop = 100,col = "white")
map('state', fill = FALSE, xlim = c(0, 100), ylim = c(0, 100), xlab = "lon", ylab = "lat", add =T)

#map("usa", xlim=c(-126.2,-65.5), ylim=c(23.8,50.6),add=T)
#map("state", xlim=c(-126.2,-65.5), ylim=c(23.8,50.6),add=T, boundary = F, interior = T, lty=2)
map("state", region="california", fill=T, add=T)

library(maps)

map("state", interior = FALSE)

map("state", boundary = FALSE, col="gray", add = TRUE)
map("state", region="california", fill=T, add=T,col='lightgray')

points(-118.5, 35.7,  pch = 21,bg='orange',col='orange')
