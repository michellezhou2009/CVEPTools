#' Construct the complete data including the observed entries and missing entries within a specifed window for regular varieties
#'
#' @param dat a data frame with the following variable names: Variety, Year, DM_Trait
#' @param w window size for the missing entries
#' @param year.list a vector of year numbers spanned over in the CVEP data
#' @param type.prev a character, "year" or "history", to indicate what type of previous performance should be evaluated
#'
#' @return a data frame which includes Variety, Year, VarietybyYear, DM_Trait, delta, New, DMtraitPrev
#'
#' @import dplyr
#'
#' @export
#'
#' @examples
compRdata_fun <- function(dat,w=Inf,year.list,type.prev="year"){

  if (any(!is.element(c("Variety","Year","DM_Trait"),colnames(data))))
    stop("Data should have the variable names: Variety, Year, DM_Trait.")

   if (!type.prev %in% c("year","history"))
    stop("type.prev is either \"year\" or \"history\".")
  ## define variable `New`
  dat = dat%>%arrange(Year)
  dat$New <- 1*(!duplicated(dat$Variety))
  dat$delta <- 1
  YearAll <- dat$Year # observed year
  if (is.infinite(w)){
    YearFull = year.list # if w=Inf, we include all the years after entry
  } else {
    YearL <- unlist(sapply(YearAll, function(x){seq(x,x+w-1,1)}, simplify = FALSE))
    YearFull = sort(unique(YearL))
  }
  YearFull = YearFull[which(YearFull %in% year.list)] # delete the years after the end of the study
  compdata = data.frame(
    "Variety" = unique(dat$Variety),
    "Year" = YearFull
  ,stringsAsFactors = F) %>%
    mutate(VarietybyYear=paste(Variety,Year,sep=":")) %>%
    left_join(dat)
  compdata$delta[is.na(compdata$delta)]=0
  compdata$New[is.na(compdata$New)]=0
  # if (any(YearFull[-1] - YearFull[-length(YearFull)]!=1)) compdata$New[which(YearFull[-1] - YearFull[-length(YearFull)]!=1)+1] = 1
  compdata = data.frame(compdata,
                   "DMtraitPrev" = sapply(1:length(YearFull),function(k){
                     if (k==1) return(NA) else {
                       tmpdat = compdata[1:(k-1),] %>% filter(!is.na(DM_Trait))
                       return(ifelse(type.prev=="year",
                                     tmpdat$DM_Trait[nrow(tmpdat)],
                                     mean(tmpdat$DM_Trait)))
                     }
                   })

  )
  compdata
}

#' Construct the complete data including the observed entries and missing entries within a specifed window for control varieties
#'
#' @param dat a data frame with the following variable names: Variety, Cycle, DM_Trait
#' @param w window size for the missing entries
#' @param cycle.list a vector of cycle numbers spanned over in the CVEP data
#' @param type.prev a character, "cycle" or "history", to indicate what type of previous performance should be evaluated
#'
#' @return a data frame which includes Variety, Cycle, VarietybyCycle, DM_Trait, delta, New, DMtraitPrev
#'
#' @import dplyr
#'
#' @export
#'
#' @examples
compCdata_fun <- function(dat,w=Inf,cycle.list,type.prev="cycle"){

  if (any(!is.element(c("Variety","Cycle","DM_Trait"),colnames(data))))
    stop("Data should have the variable names: Variety, Cycle, DM_Trait.")
  if (!type.prev %in% c("cycle","history"))
    stop("type.prev is either \"cycle\" or \"history\".")

  ## define variable `New`
  dat = dat%>%arrange(Cycle)
  dat$New <- 1*(!duplicated(dat$Variety))
  dat$delta <- 1
  CycleAll <- dat$Cycle # observed year
  if (is.infinite(w)){
    CycleFull = cycle.list # if w=Inf, we include all the years after entry
  } else {
    CycleL <- unlist(sapply(CycleAll, function(x){seq(x,x+w-1,1)}, simplify = FALSE))
    CycleFull = sort(unique(CycleL))
  }
  CycleFull = CycleFull[which(CycleFull %in% cycle.list)] # delete the years after the end of the study
  compdata = data.frame(
    "Variety" = unique(dat$Variety),
    "Cycle" = CycleFull
    ,stringsAsFactors = F) %>%
    mutate(VarietybyCycle=paste(Variety,Cycle,sep=":")) %>%
    left_join(dat)
  compdata$delta[is.na(compdata$delta)]=0
  compdata$New[is.na(compdata$New)]=0
  # if (any(CycleFull[-1] - CycleFull[-length(CycleFull)]!=1)) compdata$New[which(CycleFull[-1] - CycleFull[-length(CycleFull)]!=1)+1] = 1
  compdata = data.frame(compdata,
                        "DMtraitPrev" = sapply(1:length(CycleFull),function(k){
                          if (k==1) return(NA) else {
                            tmpdat = compdata[1:(k-1),] %>% filter(!is.na(DM_Trait))
                            return(ifelse(type.prev=="year",
                                          tmpdat$DM_Trait[nrow(tmpdat)],
                                          mean(tmpdat$DM_Trait)))
                          }
                        })

  )
  compdata
}

#' Summarize the results from an object form \code{lmer}
#'
#' @param fit.lmer a \code{\link{lmer.object}}
#' @param effects.list lists of year, location and variety
#' @param yes.IPW a logical varible indicating whether the model is fitted with weights
#'
#' @return a list including the estimates of fixed effects (includes the year effects and location effects) and random effects (including the variety, variety by year, and variety by location effects) as well as the estimates of the variance parameter (and standard deviation parameter) of random effects.
#'
#' @export
#'
#' @examples
summary.fit <- function(fit.lmer,effects.list,yes.IPW=T){

  fixeff.est = fixef(fit.lmer)
  mu.est = fixeff.est["(Intercept)"]
  year.est = fixeff.est[grep("Year",names(fixeff.est))]
  cy = matrix(rep(1,length(year.est)),ncol=1)
  year.est = c(year.est,-sum(year.est))
  loc.est = fixeff.est[grep("Loc",names(fixeff.est))]
  cl = matrix(rep(1,length(loc.est)),ncol=1)
  loc.est = c(loc.est,-sum(loc.est))
  fixeff.cov = vcov(fit.lmer)
  year.cov = fixeff.cov[grep("Year",names(fixeff.est)),grep("Year",names(fixeff.est))]
  loc.cov = fixeff.cov[grep("Loc",names(fixeff.est)),grep("Loc",names(fixeff.est))]
  mu.se = sqrt(fixeff.cov["(Intercept)","(Intercept)"])
  year.se = c(sqrt(diag(year.cov)),sqrt(as.vector(t(cy)%*%year.cov%*%cy)))
  loc.se = c(sqrt(diag(loc.cov)),sqrt(as.vector(t(cl)%*%loc.cov%*%cl)))

  out.fixeff = data.frame(
    "est"=c(mu.est,year.est,loc.est),
    "se" = c(mu.se,year.se,loc.se)
      ) %>% mutate("lwr"=est-1.96*se,"upr"=est+1.96*se)
  rownames(out.fixeff) = c("Mean",paste0("Year",effects.list$Year),paste0("Loc",effects.list$Loc))
  write.csv(out.fixeff,paste0("YearLoc_effects_",ifelse(yes.IPW,"IPW","unwgt"),".csv"))

  out = VarCorr(fit.lmer) %>% as.data.frame()
  var.raneff = out[match(c("Variety","Variety:Year","Variety:Loc","Residual"),out$grp),c("vcov","sdcor")]
  rownames(var.raneff) = c("Variety","Variety:Year","Variety:Loc","Residual")
  write.csv(var.raneff,paste0("RanEff_Var_",ifelse(yes.IPW,"IPW","unwgt"),".csv"))

  raneff.est = ranef(fit.lmer)
  ## Variety by Year Effects
  raneff.vy = raneff.est$`Variety:Year`
  raneff.vy = cbind(raneff.vy,do.call(rbind,lapply(rownames(raneff.vy),function(u){unlist(strsplit(u,split=":"))})))
  colnames(raneff.vy) = c("Est","Variety","Year")
  raneff.vy = split(raneff.vy,raneff.vy$Variety)
  out.vy = do.call(rbind,lapply(raneff.vy,function(x){
    out = rep(NA,length(effects.list$Year));names(out) = effects.list$Year
    out[match(x$Year,names(out))] = x$Est
    out
  }))
  ## Variety by Loc Effects
  raneff.vl = raneff.est$`Variety:Loc`
  raneff.vl = cbind(raneff.vl,do.call(rbind,lapply(rownames(raneff.vl),function(u){unlist(strsplit(u,split=":"))})))
  colnames(raneff.vl) = c("Est","Variety","Loc")
  raneff.vl = split(raneff.vl,raneff.vl$Variety)
  out.vl = do.call(rbind,lapply(raneff.vl,function(x){
    out = rep(NA,length(effects.list$Loc));names(out) = effects.list$Loc
    out[match(x$Loc,names(out))] = x$Est
    out
  }))
  out.raneff = data.frame("Variety" = raneff.est$Variety[match(effects.list$Variety,rownames(raneff.est$Variety)),1])
  out.raneff = cbind(out.raneff,out.vy[match(effects.list$Variety,rownames(out.vy)),],out.vl[match(effects.list$Variety,rownames(out.vl)),])
  write.csv(out.raneff,paste0("RanEff_estimates_",ifelse(yes.IPW,"IPW","unwgt"),".csv"))

  return(list("fixeff"=out.fixeff,"raneff"=out.raneff,"raneff.var"=var.raneff))
}
