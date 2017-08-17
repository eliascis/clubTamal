#' @title Transform Panel Data
#' @description panel.transform converts the a panel data frame to time-demeand data or to first differenced data.
#' @param formula a formula, indicating the covariates.
#' @param data a data frame
#' @param model a model, "fd" or "within"
#' @param index id and group structure, e.g. c("id","year")
#' @details
#' convert data frame into transformed data
#' @return data.frame
#' @author Elías Cisneros <ec@elias-cisneros.de>
#' @importFrom stats lm
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom utils head
#' @import plm
#' @export


panel.transform<-function(
  formula,
  data,
  model="within",
  index=c("id","year")
){
  # formula<-y~x+factor(year)
  # data<-spd4testing()
  # model="within"
  # index=c("id","year")

  ##setup
  f<-formula
  pf<-pFormula(f)
  data<-pdata.frame(data,index)

  ##restricitons
  if ((model %in% c("fd","within"))!=T){
    stop('panel.transform transforms models of c("fd","within")')
  }

  ##transformed data
  #index

  #reduced data
  # x<-d[,c(pindex,as.character(formula[2]))]
  x<-data[,index]
  x$rowid<-rownames(x)
  rd<-x


  #model frame
  x<-model.frame(f,data)
  mf<-x
  #

  #model response
  x<-pmodel.response(object=pf,data=data,model)
  x<-data.frame(x)
  names(x)<-names(mf)[1]
  mr<-x

  #model matrix - transformed data of explanatory variables
  x<-model.matrix(pf,data,model)
  # n<-names(e$aliased[e$aliased==F])
  # i<-match(n,colnames(x))
  # i<-i[!is.na(i)] #avoid omited constant variables in new data set
  # x<-x[,i]
  # n<-colnames(x)
  # x<-data.frame(x)
  # names(x)<-n
  mx<-x

  #transformed data = model response + model matrix
  x<-cbind(mr,mx)
  # i<-match("rowid",names(x))
  # i<-c(i,match("(intercept)",names(x)))
  # i<-i[!is.na(i)]
  # if (length(i)>0){
  #   x<-x[,-i]
  # }
  # coefnames<-names(x)[-1]
  # x<-data.frame(x) #converts names to be usable in regression
  mt<-x

  return(mt)

}
