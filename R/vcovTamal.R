#' @title Group clustering of panel data esitmates
#' @description vcovTamal, creates a robust covariance variance matrix clustered at group level for \code{lm} objects.
#' Clustering is based on Angist and Pischke (2009) resulting in Stata(TM) like standard errors of
#' regression coefficients.
#' @param estimate an object of class \code{"plm"} estimated with by one of the methods \code{within}, \code{fd}.
#' @param data a data.frame object, that must be the data.frame used to create \code{estimate} object.
#' Can be a data.frame or a pdata.frame object.
#' @param groupvar a string indicating a column in \code{data} to indexes the group structure.
#' @param byhand logical, if TRUE, the clustered covariance matrix is calculated by formulas
#' without using the \code{multiwayvcov} package. This method is outdated.
#' @details
#' See package vignette.
#' @return A matrix containing the covariance matrix estimate
#' @author Elías Cisneros <ec@elias-cisneros.de>
#' @example man/eg.vcovTamal.R
#' @references
#' Angrist, J. D. & Pischke, J.-S. Mostly harmless econometrics: An empiricist's companion Princeton university press, 2009
#' Arai, M. Cluster-robust standard errors using R, 2015.
#' See also: https://thetarzan.wordpress.com/2011/06/11/clustered-standard-errors-in-r/.
#' @importFrom stats lm
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom utils head
#' @importFrom lmtest coeftest
#' @import plm
#' @import sandwich
#' @import multiwayvcov
#' @note The function uses the \code{estimate} and \code{data} to construct a transfomed data frame
#' - first differenced or time-demeaned. The esitmate is then re-estimated by OLS with \code{lm} in order to extract the
#' covarinace matrix and correct for degrees of freedom and use a sandwich estimate with the group structure.
#' @export


vcovTamal<-function(
  estimate,
  data,
  groupvar,
  byhand=F
){

  #example
  # data<-spd4testing(missingX=T)
  # data<-pdata.frame(data)
  # f<-formula(y~x+factor(year))
  # estimate<-plm(f,data,method="within")

  #setup
  model<- estimate$args$model

  ##restricitons
  if (model %in% c("within")){
    stop('better use the very fast and powerful lfe::felm')
  }
  if ((model %in% c("fd","within"))!=T){
    stop('vcovTamal only applys to plm models of c("fd")')
  }

  ##convert plm to lm
  et<-plm2lm(estimate=estimate,reestimate=T)

  ##group index
  if (is.numeric(groupvar)==T){
    stop("groupvar is numeric. not yet implmeted")
  }
  if (is.character(groupvar)==T & is.data.frame(data)==T){
    x<-data
    x$rowid<-rownames(x)
    x<-x[x$rowid %in% rownames(et$model),c("rowid",groupvar)]
    x<-x[,groupvar]
    x<-as.integer(x)
    gx<-x
  }

  ##cluster with multiwayvcov
  if (model=="within"){
    vcovCL<-cluster.vcov(et,gx,stata_fe_model_rank=T)
  }
  if (model=="fd"){
    vcovCL<-cluster.vcov(et,gx,stata_fe_model_rank=F)
  }

  ##cluster by hand (outdated)
  if (byhand==T){
    ##cluster
    #weighting of covaraince matrix
    M<-length(unique(gx))
    N<-length(gx)
    K<-et$rank
    if (model=="fd"){
      dfcw <- 1
    }
    if (model=="within"){
      dfcw <- et$df / (et$df - (M-1))  # (N-K)/( (N-K)-(M-1) )
    }

    #degree of freedom correction
    dfc<-(M/(M-1))*((N-1)/(N-K)) # (M/(M-1)) * ((N-1)/(N-K))

    #group residuals uj = Xj * ej
    u<-apply(
      estfun(et),
      2,
      function(x) {
        tapply(x, gx, sum)
      }
    )
    #sandwich estimator
    # meat<-crossprod(u)/N
    # sx <- summary.lm(et)
    # bread<-sx$cov.unscaled * as.vector(sum(sx$df[1:2]))
    # vcovCL<- dfc*1/N*bread%*%meat%*%bread
    vcovCL<- dfc * sandwich(x=et, meat.=crossprod(u)/N) * dfcw
  }

  ##test
  # coeftest(et, vcovCL)

  ##out
  return(vcovCL)
}




# https://rdrr.io/rforge/plm/src/R/plm.R
# https://r-forge.r-project.org/scm/viewvc.php/pkg/R/pfunctions.R?view=markup&root=plm&sortdir=down

