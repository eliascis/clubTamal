#' @title Transform plm to lm object
#' @description plm2lm, converts the a \code{plm} object to an \code{lm} object.
#' The information saved in a \code{plm} object is used to obtain an transformed dataset
#' (e.g. time-demeaned or first differences)
#' in order to restimate an OLS model with \code{lm}.
#' @param estimate an object of class \code{"plm"}.
#' @param reestimate logical, if TRUE the panel data model is reestimated with a linear regression (\code{lm}).
#' @details
#' estimate - The \code{plm} object estimated by one of the methods \code{within}, \code{fd}.
#' reestimate - If FALSE, an error displays - A conversion without reestimation is currently not implmented.
#' @return A matrix containing the covariance matrix estimate
#' @author Elías Cisneros <ec@elias-cisneros.de>
#' @example man/eg.vcovTamal.R
#' @importFrom stats lm
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom utils head
#' @import plm
#' @note Reestimation of results leads to an increase in processing requairemets.
#' @export


plm2lm<-function(estimate,reestimate=T){

  ###conversion by reestimation
  #example
  # d<-spd4testing(missingX=T)
  # f<-formula(y~x+factor(year))
  # e<-plm(f,d,model="within")

  ##setup
  e<-estimate
  formula<-e$formula
  model<-e$args$model

  ##restricitons
  if ((model %in% c("fd","within"))!=T){
    stop('vcovTamal only applys to plm models of c("fd","within")')
  }


  ##transformed data
  #index
  pindex<-names(attributes(e$model)$index)

  #reduced data
  # x<-d[,c(pindex,as.character(formula[2]))]
  x<-attributes(e$model)$index
  x$rowid<-as.integer(rownames(x))
  rd<-x

  #model frame
  x<-e$model
  mf<-x

  #model response
  x<-pmodel.response(formula,mf,model)
  x<-data.frame(x)
  names(x)<-names(mf)[1]
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
  mx<-x

  #transformed data = model response + model matrix
  x<-cbind(mr,mx)
  i<-match("rowid",names(x))
  i<-c(i,match("(intercept)",names(x)))
  i<-i[!is.na(i)]
  if (length(i)>0){
    x<-x[,-i]
  }
  coefnames<-names(x)[-1]
  x<-data.frame(x) #converts names to be usable in regression
  mt<-x


  ##new estimation on transformed data
  #formula
  if (model=="fd"){
    ft<-paste(names(mt)[1], "~ ", paste(names(mt)[-1],collapse="+"))
  }
  if (model=="within") {
    ft<-paste(names(mt)[1], "~ 0+", paste(names(mt)[-1],collapse="+"))
  }
  #estimation
  et<-lm(ft,mt)
  #estimation names correction
  if (model=="fd"){
    names(et$coefficients)<-c("(intercept)",coefnames)
  }
  if (model=="within"){
    names(et$coefficients)<-coefnames
  }


  ###conversion without reestimation
  if (reestimate!=T){
    stop("A simple conversion of plm to lm is still pending")

    z<-list()
    z$coefficients<-e$coefficients
    # z$na.action <- NULL
    z$residuals <- e$residuals
    z$effects<-NA
    z$rank <- length(e$coefficients)
    z$fitted.values <- e$model[,1]-e$residuals
    z$assign<-e$assign
    # z$qr<-list(
    #   qr = model.matrix(
    #     e$formula,
    #     mf,
    #     model
    #   )
    #   pmodel.response(e$model)
    z$qr<-NA
    z$df.resiidual<-e$df.residual
    z$contrast<-e$contrasts
    z$xlevels<-NA
    z$call<-e$call
    z$terms<-NA
    z$model<-e$model
    class(z)<-"lm"
    names(z)

    summary(z)
    et<-z
  }

  ###out
  return(et)
}

