#'Interpolating points to data of list
#'@param collect data to be interpolated
#'@importFrom myfs debugText
#'@ export
#'
all_interpolate<-function(collect){
  incollect<-collect
  for (l in 1:length(collect)) {
    inter.oricord<-voronoi_interpo(collect[[l]][[2]], 15)
    incollect[[l]][[2]]<-rbind(incollect[[l]][[2]], inter.oricord)
    incollect[[l]][[1]]<-nrow(incollect[[l]][[2]])
    myfs::debugText(l, incollect[[l]][[1]])
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

      vics<-get_vicinity(dist, i, nvics)

      element[vics]<-element[vics]+1

      #vics.oricord<-voronoiProcess(vics.line, figure)
      vics.oricord<-voronoiBorder(vics, figure)[[1]]

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

#'Add points at the vertexes of Voronoi region
#'@param vics center and neighbor points
#'@param figure data set to be interpolated
#'@importFrom deldir deldir
#'@return added points
#'
voronoi_border<-function(vics, figure){


  vics_pca<-prcomp(figure[vics,])

  res<-deldir::deldir(vics_pca$x[,1], vics_pca$x[,2])

  tiles<-tile.list(res)

  insecs<-cbind(tiles[[1]][["x"]], tiles[[1]][["y"]])

  exist<-exist_convexhull_check(vics_pca, insecs)


  vics.oricord<-originCoodinate(vics.pca, insecs[which(exist==T), ])

  return(list(oricord=vics.oricord, pca.inter=insecs[which(exist==T), ]))

}



#'Distinguishing whether passed points are in convex.
#'@param rpca points consisting a convex
#'@param insecs points to be distinguished
#'@importFrom grDevices chull
#'
#'
#'
exist_convexhull_check<-function(rpca, insecs){

  chul<-grDevices::chull(rpca[["x"]][,1:2])

  exist<-sapply(1:nrow(insecs), function(i){
    sides<-chul[which(rpca[["x"]][chul,1]>=insecs[i,1])] %>%
      sapply(., function(k)convex_hull_vertx(chul, k))
    if(length(sides)==0){return(F)}
    else{cross.side<-sidesSet(sides)}
    hline<-matrix(c(insecs[i,], max(rpca[["x"]][chul,1][which(rpca[["x"]][chul,1]>=insecs[i,1])]), insecs[i,2]), 2, 2, byrow=T)
    #debugText(hline)
    c.ncross<-convex_hull_check(rpca, hline, t(cross.side))
    #debugText(c.ncross)
    if(length(which(c.ncross==T)) %% 2 != 0){return(T)}
    else{return(F)}

  })

  return(exist)

}
