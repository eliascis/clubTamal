##packages
library(foreign)
library(clubTamal)
library(lmtest)
library(plm)
library(clubSandwich)
library(spd4testing)

##data
d<-spd4testing()
d

#manipulate data
# d<-pdata.frame(d)

##formula
f<- formula(y ~ x + factor(year))

##standard estimation
e<-plm(formula=f,data=d,model="fd")
summary(e)
e<-plm(formula=f,data=d,model="within")
summary(e)

##clustering
#no clustering
v<-e$vcov
coeftest(e)
#clustering at id level with plm-package
v<-vcovHC(e,type="HC1",cluster="group",tol=1*10^-20 )
coeftest(e,v)

##clustering at group level with clubSandwich-package
v<-vcovCR(e,d$gid,type="CR1")
coef_test(e,vcov=v,test="naive-t")
coeftest(e,v)

##clustering at group level with clubTamal
v<-vcovTamal(estimate=e,data=d,groupvar="gid")
coeftest(e,v)


##clustering at group level with STATA (TM)
#export data
tdir<-tempdir()
write.csv(d,file.path(tdir,"spd.csv"))
#write stata file
sink(file.path(tdir,"clusterStata.do"))
cat(paste0('insheet using "',file.path("spd.csv"),'", clear'))
cat("\n")
cat('xtset id year')
cat("\n")
cat('xtreg y x i.year, fe')
cat("\n")
cat('xtreg y x i.year, fe vce(cluster gid)')
cat("\n")
cat('reg d.y d.x i.year')
cat("\n")
cat('reg d.y d.x i.year, vce(cluster gid)')
sink()
## RUN Stata file

