#' @title Group clustering of panel data esitmates
#' @description vcovTamal, creates a robust covariance variance matrix clustered at group level for \code{lm} objects.
#' Clustering is based on Angist and Pischke (2009) resulting in Stata(TM) like standard errors of
#' regression coefficients.
#' @param estimate an object of class \code{"plm"} estimated with by one of the methods \code{within}, \code{fd}.
#' @param data the data frame object used to create \code{estimate} object. Can be a data.frame or a pdata.frame object.
#' @param groupvar a string indicating a column in \code{data} to indexes the group structure.
#' @details
#' Clustering \cite{plm} objects estimated by a "within" model: \cr
#' Details forthcomming... \cr
#' @details
#' Clustering \cite{plm} objects estimated by a "fd" model: \cr
#' Details forthcomming... \cr
#' @return A matrix containing the covariance matrix estimate
#' @author Elías Cisneros <ec@elias-cisneros.de>
#' @example man/eg.vcovTamal.R
#' @references Angrist, J. D. & Pischke, J.-S. Mostly harmless econometrics: An empiricist's companion Princeton university press, 2009
#' @importFrom stats lm
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom utils head
#' @importFrom lmtest coeftest
#' @import plm
#' @import sandwich
#' @note The function uses the \code{estimate} and \code{data} to construct a transfomed data frame
#' - first differenced or time-demeaned. The esitmate is then re-estimated by OLS with \code{lm} in order to extract the
#' covarinace matrix and correct for degrees of freedom and use a sandwich estimate with the group structure.
#' @export


vcovTamal<-function(
  estimate,
  data,
  groupvar
){

  # estimate<-d
  # data<-d
  # groupvar<-"gid"

  #setup
  e<-estimate
  d<-data
  vc<-e$vcov
  formula<-e$formula
  model<-e$args$model

  #restricitons
  if ((model %in% c("fd","within"))!=T){
    stop('vcovTamal only applys to plm models of c("fd","within")')
  }

  #index
  pindex<-names(attributes(e$model)$index)

  #reduced data
  # x<-d[,c(pindex,as.character(formula[2]))]
  x<-d[,c(pindex,groupvar)]
  rownames(x)<-1:nrow(x)
  x$rowid<-as.integer(rownames(x))
  head(x)
  rd<-x

  #model frame
  if (class(d)[1]=="pdata.frame"){
    d.pdata.frame<-d
  }
  if (class(d)[1]=="data.frame"){
    d.pdata.frame<-pdata.frame(d)
  }
  x<-model.frame(formula,data=d.pdata.frame)
  # x$rowid<-as.integer(rownames(x))
  # head(x)
  mf<-x

  #model response
  x<-pmodel.response(formula,mf,model)
  x<-data.frame(x)
  names(x)<-names(mf)[1]
  x$rowid<-as.integer(rownames(x))
  mr<-x

  #model matrix - transformed data of explanatory variables
  x<-model.matrix(formula,mf,model)
  n<-names(e$aliased[e$aliased==F])
  i<-match(n,colnames(x))
  i<-i[!is.na(i)] #avoid omited constant variables in new data set
  x<-x[,i]
  n<-colnames(x)
  x<-data.frame(x)
  names(x)<-n
  x$rowid<-as.integer(rownames(x))
  mx<-x

  #transformed data
  x<-merge(mr,mx,by="rowid",all.x=T)
  i<-match("rowid",names(x))
  i<-c(i,match("(intercept)",names(x)))
  i<-i[!is.na(i)]
  if (length(i)>0){
    x<-x[,-i]
  }
  coefnames<-names(x)[-1]
  x<-data.frame(x) #confert names to be usable in regression
  mt<-x

  ##new estimation on transformed data
  if (model=="fd"){
    ft<-paste(names(mt)[1], "~ ", paste(names(mt)[-1],collapse="+"))
  }
  if (model=="within") {
    ft<-paste(names(mt)[1], "~ 0+", paste(names(mt)[-1],collapse="+"))
  }
  et<-lm(ft,mt)
  if (model=="fd"){
    names(et$coefficients)<-c("(intercept)",coefnames)
  }
  if (model=="within"){
    names(et$coefficients)<-coefnames
  }

  #group index
  x<-rd
  x<-x[x$rowid %in% mx$rowid,c("rowid",pindex,groupvar)]
  head(x)
  nrow(x)
  x<-x[,groupvar]
  gx<-x

  ##cluster
  #degrees of freedom (robust)
  M<-length(unique(gx))
  N<-length(gx)
  K<-et$rank
  if (model=="fd"){
    dfcw <- 1
  }
  if (model=="within"){
    dfcw <- et$df / (et$df - (M-1))  # dfcw<- (N-K)/( (N-K)-(M-1) )
  }

  clx <- function(fm=et, dfcw=dfcw, cluster=gx){
    # fm<-et
    # dfcw<-dfcw
    # cluster<-gx

    M<-length(unique(cluster))
    N<-length(cluster)
    dfc<-(M/(M-1))*((N-1)/(N-fm$rank)) # dfc <- (M/(M-1)) * ((N-1)/(N-K))
    u<-apply(
      estfun(fm),
      2,
      function(x) {
        tapply(x, cluster, sum)
      }
    )
    # meat<-crossprod(u)/N
    # sx <- summary.lm(fm)
    # bread<-sx$cov.unscaled * as.vector(sum(sx$df[1:2]))
    # vcovCL<- dfc*1/N*bread%*%meat%*%bread
    vcovCL<- dfc * sandwich(x=fm, meat.=crossprod(u)/N) * dfcw
    # coeftest(fm, vcovCL)
    return(vcovCL)
  }
  vcn<-clx(et, dfcw, gx)

  ##out
  return(vcn)
}




# https://rdrr.io/rforge/plm/src/R/plm.R
# https://r-forge.r-project.org/scm/viewvc.php/pkg/R/pfunctions.R?view=markup&root=plm&sortdir=down

