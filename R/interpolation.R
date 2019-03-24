#'
#'
#'
all_interpolate<-function(collect){
  incollect<-collect
  for (l in 1:length(collect)) {
    inter.oricord<-voronoiInterpo(collect[[l]][["noizyX"]], 15)
    incollect[[l]][["noizyX"]]<-rbind(incollect[[l]][["noizyX"]], inter.oricord)
    incollect[[l]][["nsample"]]<-nrow(incollect[[l]][["noizyX"]])
    debugText(l, incollect[[l]][["nsample"]])
  }

  return(incollect)
}
