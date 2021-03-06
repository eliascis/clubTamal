% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcovTamal.R
\name{vcovTamal}
\alias{vcovTamal}
\title{Group clustering of panel data esitmates}
\usage{
vcovTamal(estimate, data, groupvar, byhand = F)
}
\arguments{
\item{estimate}{an object of class \code{plm} estimated with the \code{fd} option in methods.}

\item{data}{a data.frame object, that must be the data.frame used to create the \code{plm} object.
Can be a data.frame or a pdata.frame object.}

\item{groupvar}{a string indicating a column in \code{data} that indicates the group structure.}

\item{byhand}{logical, if TRUE, the clustered covariance matrix is calculated by formulas
without using the \code{multiwayvcov} package. This method is outdated.}
}
\value{
A matrix containing the covariance matrix estimate
}
\description{
vcovTamal, creates a robust covariance variance matrix clustered at group level for first difference estimations of the \code{lm::plm} command.
Clustering is based on Angist and Pischke (2009) resulting in Stata(TM) like standard errors of
regression coefficients.
}
\details{
See package vignette.
}
\note{
The function uses the \code{plm} output and the \code{data} to construct a transfomed (first differenced) data frame.
The esitmate is then re-estimated by OLS with \code{lm} in order to extract the
covarinace matrix and correct for degrees of freedom and use a sandwich estimate with the group structure.
The implementation for Fixed Effects (within) estimates is outdated. Ihe command \code{felm} of the package \code{lfe} is able construct directly estimates with multiple cluster with multiple fixed effects.
}
\examples{
##packages
library(foreign)
library(clubTamal)
library(lmtest)
library(plm)
library(clubSandwich)
library(spd4testing)
library(multiwayvcov)

##data
d<-spd4testing()
d

##formula
f<- formula(y ~ x + factor(year))

##standard estimation
e<-plm(formula=f,data=d,model="fd")
summary(e)
\dontrun{
  # e<-plm(formula=f,data=d,model="within")
  # summary(e)
}

##clustering
#no clustering
v<-e$vcov
coeftest(e)
#clustering at id level with plm-package
v<-vcovHC(e,type="HC1",cluster="group",tol=1*10^-20 )
coeftest(e,v)

##clustering at group level with clubSandwich-package
##only for Fixed Effects regressions
v<-vcovCR(e,d$gid,type="CR1")
coef_test(e,vcov=v,test="naive-t")
coeftest(e,v)

##clustering at group level with clubTamal
v<-vcovTamal(estimate=e,data=d,groupvar="gid")
coeftest(e,v)

##clustering with clubTamal but without using the multivcov package
\dontrun{
v<-vcovTamal(estimate=e,data=d,groupvar="gid",byhand=T)
coeftest(e,v)
}

##clustering at group level with STATA (TM)
#export data
tdir<-tempdir()
write.csv(d,file.path(tdir,"spd.csv"))
#write stata file
sink(file.path(tdir,"clusterStata.do"))
cat(paste0('insheet using "',file.path("spd.csv"),'", clear'))
cat("\\n")
cat('xtset id year')
cat("\\n")
cat('xtreg y x i.year, fe')
cat("\\n")
cat('xtreg y x i.year, fe vce(cluster gid)')
cat("\\n")
cat('reg d.y d.x i.year')
cat("\\n")
cat('reg d.y d.x i.year, vce(cluster gid)')
sink()
## RUN Stata file


}
\references{
Angrist, J. D. & Pischke, J.-S. Mostly harmless econometrics: An empiricist's companion Princeton university press, 2009
Arai, M. Cluster-robust standard errors using R, 2015.
See also: https://thetarzan.wordpress.com/2011/06/11/clustered-standard-errors-in-r/.
}
\author{
Elías Cisneros <ec@elias-cisneros.de>
}
