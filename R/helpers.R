


compRdata_fun <- function(dat,w=Inf){
  ## define variable `New`
  dat$New <- 1*(!duplicated(dat$Variety))
  YearAll <- dat$Year # observed year
  if (is.infinite(w)){
    YearFull = year.list # if w=Inf, we include all the years after entry
  } else {
    YearL <- unlist(sapply(YearAll, function(x){seq(x,x+w-1,1)}, simplify = FALSE))
    YearFull = sort(unique(YearL))
  }
  YearFull = YearFull[which(YearFull %in% year.list)] # delete the years after the end of the study
  compdata = data.frame(
    "Variety"=unique(dat$Variety),
    "Year" = YearFull,
    "Cycle"=cycle_year%>%filter(Year%in%YearFull)%>%select(Cycle)
  ) %>% 
    dplyr::mutate(VarietybyYear=paste(Variety,Year,sep=":")) %>% 
    dplyr::left_join(dat)
  compdata$delta[is.na(compdata$delta)]=0
  compdata$New[is.na(compdata$New)]=0
  if (any(YearFull[-1] - YearFull[-length(YearFull)]!=1)) compdata$New[which(YearFull[-1] - YearFull[-length(YearFull)]!=1)+1] = 1 
  compdata = cbind(compdata,
                   do.call(rbind,
                           lapply(1:length(YearFull),function(k){
                             if (k==1) { 
                               data.frame("DMtrait_Prev"=NA,
                                          "DMtrait_PrevAll"=NA,
                                          "NumYear"=0) 
                             } else {
                               tmpdat = compdata[1:(k-1),] %>% dplyr::filter(!is.na(DM_Trait))
                               data.frame("DMtrait_Prev" = tmpdat$DM_Trait[nrow(tmpdat)],
                                          "DMtrait_PrevAll" = mean(tmpdat$DM_Trait),
                                          "NumYear" = sum(compdata$delta[1:(k-1)])
                               )
                             } # end for else 
                           })
                   )
  )
  compdata
}
