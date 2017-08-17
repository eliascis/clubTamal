d<-spd4testing()
d

d.t<-panel.transform(y~x+factor(year),d,model="within")
