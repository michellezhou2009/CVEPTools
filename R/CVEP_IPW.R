#' Fit a linear mixed effects model with inverse probability weighting
#'
#' @param data a data frame with the following variable names: Variety, Year, Loc, Rep, Cycle, Checks, TT_Trait,t DM_Trait
#' @param control a list of control parameters including \code{wyear} which specifies the window size in years of the missing entries for regular varieties and \code{wcycle} which specifies the window size in cycles of the missing entries for control varieties
#' @param yes.unwgt a logical variable indicating whether to fit an unweighted linear mixed effects model
#'
#' @return a list of data used to fit the linear mixed effects model, the fitted model and summarized results
#'
#' @import dplyr
#' @import lme4
#'
#' @export
#'
#' @examples
CVEP_IPW <- function(data,control=list(wyear=3,wcycle=3),yes.unwgt=F){

  if (any(!is.element(c("Variety","Year","Loc","Rep","Cycle","Checks", "TT_Trait","DM_Trait"),colnames(data))))
    stop("Data should have the variable names: Variety, Year, Loc, Rep, Cycle, Checks, TT_Trait,t DM_Trait.")
  ## ??? need to check the characters of "Check"

  data = data %>% mutate(
    Variety = as.character(data$Variety)
  )
  year.list = sort(unique(data$Year))
  cycle.list = sort(unique(data$Cycle))
  cycle_year = data %>% group_by(Year) %>% summarise(Cycle=unique(Cycle))
  variety.list = unique(data$Variety)
  variety.control.list = unique(data$Variety[data$Checks=="Yes"])
  variety.regular.list = setdiff(variety.list,variety.control.list)
  loc.list = unique(data$Loc)

  data_R = data %>% filter(Variety %in% variety.regular.list) %>%
    mutate(VarietybyYear=paste(Variety,Year,sep=":"))
  data_C = data %>% filter(Variety %in% variety.control.list) %>%
    mutate(VarietybyCycle=paste(Variety,Cycle,sep=":"))

  ### ==== compute the probability whether each regular variety is tested at the years following entering the program
  DMdat = data_R %>%
    group_by(VarietybyYear) %>%
    summarise(
      Variety=as.character(unique(Variety)),
      Year=unique(Year),
      DM_Trait=unique(DM_Trait)
      )
  DMdat_Variety = split(DMdat, as.factor(DMdat$Variety))
  DMcompdata_all = lapply(DMdat_Variety,compRdata_fun,
                          w=control$wyear, year.list=year.list, type.prev="year")
  DMcompdata = do.call(rbind,DMcompdata_all)
  DMcompdata$prob_delta = DMcompdata$New
  ncase_year = DMcompdata %>% filter(New==0) %>% group_by(Year) %>% summarise(ncase=sum(delta))
  index = DMcompdata$New==0 & !DMcompdata$Year%in%ncase_year$Year[ncase_year$ncase==0]
  fit.dat = DMcompdata %>% filter(index) %>%
    mutate(Yearf=C(factor(Year),contr = sum),
           DMtraitPrev = DMtraitPrev - 1)
  fit_R = glm(delta ~ Yearf + DMtraitPrev, family=binomial, fit.dat)

  DMcompdata$prob_delta[index] = fit_R$fitted.values
  data_R = data_R %>% mutate(prob_delta=DMcompdata$prob_delta[match(VarietybyYear,DMcompdata$VarietybyYear)])

  ### ==== compute the probability whether each control variety is tested at the years following entering the program
  DMdat = data_C %>%
    group_by(VarietybyCycle) %>%
    summarise(
      Variety=as.character(unique(Variety)),
      Cycle = unique(Cycle),
      DM_Trait=mean(DM_Trait,na.rm=T)
    )
  DMdat_Variety = split(DMdat, as.factor(DMdat$Variety))
  DMcompdata_all = lapply(DMdat_Variety,compCdata_fun,
                          w=control$wcycle, cycle.list=cycle.list, type.prev="cycle")
  DMcompdata = do.call(rbind,DMcompdata_all)
  DMcompdata$prob_delta = DMcompdata$New
  ncase_cycle = DMcompdata %>% filter(New==0) %>% group_by(Cycle) %>% summarise(ncase=sum(delta))
  index = DMcompdata$New==0 & !DMcompdata$Cycle%in%ncase_cycle$Cycle[ncase_cycle$ncase==0]
  fit.dat = DMcompdata %>% filter(index) %>%
    mutate(Cyclef=C(factor(Cycle),contr = sum),
           DMtraitPrev = DMtraitPrev - 1)
  fit_C = glm(delta ~ Cyclef + DMtraitPrev, family=binomial, fit.dat)
  DMcompdata$prob_delta[index] = fit_C$fitted.values
  data_C = data_C %>% mutate(prob_delta=DMcompdata$prob_delta[match(VarietybyCycle,DMcompdata$VarietybyCycle)])

  ## Calculate the probability that the trait of a given variety at a given year, location and replicate is observed
  prob.year = data %>% group_by(Year) %>%
    summarize(prob = sum(1*!is.na(TT_Trait))/(length(unique(Variety))*length(unique(Loc))*length(unique(Rep)))
    )

  ## construct the data for fitting the linear mixed effects model with calculated weights
  newdata = rbind(
    data_R %>% select(Variety,Year,Loc,Rep,Cycle,Checks,TT_Trait, DM_Trait,prob_delta),
    data_C %>% select(Variety,Year,Loc,Rep,Cycle,Checks,TT_Trait, DM_Trait,prob_delta)
    ) %>%
    mutate(prob_year=prob.year$prob[match(Year,prob.year$Year)],
           weight = 1/(prob_delta*prob_year))

  newdata = newdata %>% filter(!is.na(TT_Trait)) %>%
    mutate(Yearf = C(factor(Year,levels=sort(unique(Year))),contr=sum),
           Locf = C(factor(Loc,levels=unique(Loc)),contr=sum))
  fit.IPW =  lmer(TT_Trait ~ Yearf + Locf + (1|Variety) + (1|Variety:Year) + (1|Variety:Loc), data=newdata, weight=weight)
  effects.list = list("Variety"=unique(as.character(newdata$Variety)),"Year"=sort(unique(newdata$Year)),"Loc"=unique(as.character(newdata$Loc)))
  if (yes.unwgt)  {
    fit.unwgt = lmer(TT_Trait ~ Year + Loc + (1|Variety) + (1|Variety:Year) + (1|Variety:Loc), data=newdata)
    return(list(fit.dat=newdata,fit.IPW=fit.IPW,summary.IPW=summary.fit(fit.IPW,effects.list), fit.unwgt=fit.unwgt,summary.unwgt=summary.fit(fit.unwgt,yes.IPW=F)))
  } else return(list(fit.dat=newdata,fit.IPW=fit.IPW,summary.IPW=summary.fit(fit.IPW,effects.list)))
}
