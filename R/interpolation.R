#'Interpolating points to data of list
#'@param collect data to be interpolated
#'@ export
#'
all_interpolate<-function(collect){
  incollect<-collect
  for (l in 1:length(collect)) {
    inter.oricord<-voronoi_interpo(collect[[l]][[2]], 15)
    incollect[[l]][[2]]<-rbind(incollect[[l]][[2]], inter.oricord)
    incollect[[l]][[1]]<-nrow(incollect[[l]][[2]])
    debugText(l, incollect[[l]][[1]])
  }

  return(incollect)
}


#'Interpolating points to one data set
#'@param figure data set
#'@param nvics the number of neighborhood points
#'
voronoiInterpo<-function(figure, nvics){

  element<-rep(0, length = nrow(figure))

  dist<-dist(figure)

  for (i in 1:nrow(figure)) {

    if(element[i]==0){

      vics<-get.vicinity(dist, i, nvics)

      vics.line<-line.vics(i, vics)

      element[vics.line]<-element[vics.line]+1

      #vics.oricord<-voronoiProcess(vics.line, figure)
      vics.oricord<-voronoiBorder(vics.line, figure)[[1]]

      if(i==1){oricord<-vics.oricord}
      else{oricord<-rbind(oricord, vics.oricord)}

    }

  }

  #debugText(element)

  return(oricord)

}

#' Choose 1-n th close points to a point
#' @param dist distant matrix
#' @param center a vector element number of center of the neighborhoods
#' @param nvic the number of neighborhood points
#' @importFrom magrittr %>%
#' @return vector element numbers the neighborhoods
#' @export
#'
get_vicinity<-function(dist, center, nvic){

  vic<-as.matrix(dist)[center, ] %>%
       order() %>%
       .[1:(nvic+1)]

  return(vic)

}
