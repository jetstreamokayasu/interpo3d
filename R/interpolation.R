#'Interpolating points to data of list
#'@param collect data to be interpolated
#' @param nvic the number of neighbor points
#'@importFrom myfs debugText
#' @return list of interpolated data sets
#'@export
#'
all_interpolate<-function(collect, nvic){
  incollect<-collect
  for (l in 1:length(collect)) {
    inter.oricord<-voronoi_interpo(collect[[l]][[2]], nvic)
    incollect[[l]][[2]]<-rbind(incollect[[l]][[2]], inter.oricord)
    incollect[[l]][[1]]<-nrow(incollect[[l]][[2]])
    cat("dataset", l, "has", incollect[[l]][[1]], "points\n")
  }

  return(incollect)
}


#'Interpolating points to one data set
#'@param figure data set
#'@param nvics the number of neighborhood points
#'@return intepolated points
#'@export
#'
voronoi_interpo<-function(figure, nvics){

  element<-rep(0, length = nrow(figure))

  dist<-dist(figure)

  for (i in 1:nrow(figure)) {

    if(element[i]==0){

      vics<-get_vicinity(dist, i, nvics)

      element[vics]<-element[vics]+1

      vics_oricord<-voronoi_border(vics, figure)[[1]]

      if(i==1){oricord<-vics_oricord}
      else{oricord<-rbind(oricord, vics_oricord)}

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
#'@importFrom deldir tile.list
#'@importFrom stats prcomp
#'@return added points
#'
voronoi_border<-function(vics, figure){


  vics_pca<-stats::prcomp(figure[vics,])

  res<-deldir::deldir(vics_pca$x[,1], vics_pca$x[,2])

  tiles<-deldir::tile.list(res)

  insecs<-cbind(tiles[[1]][["x"]], tiles[[1]][["y"]])

  exist<-exist_convexhull_check(vics_pca, insecs)

  if(length(insecs[which(exist==T), ]) > 0){

    vics_oricord<-origin_coordinate(vics_pca, insecs[which(exist==T), ], figure[vics[1],])

    return(list(oricord=vics_oricord, pca_inter=insecs[which(exist==T), ]))

  }

}



#'Distinguishing whether passed points are in convex.
#'@param rpca points consisting a convex
#'@param insecs points to be distinguished
#'@importFrom grDevices chull
#'@return True when points are in convex
#'
#'
exist_convexhull_check<-function(rpca, insecs){

  chul<-grDevices::chull(rpca[["x"]][,1:2])

  exist<-sapply(1:nrow(insecs), function(i){
    sides<-chul[which(rpca[["x"]][chul,1]>=insecs[i,1])] %>%
      sapply(., function(k)convex_hull_vertx(chul, k))
    if(length(sides)==0){return(F)}
    else{cross_side<-sides_set(sides)}
    hline<-matrix(c(insecs[i,], max(rpca[["x"]][chul,1][which(rpca[["x"]][chul,1]>=insecs[i,1])]), insecs[i,2]), 2, 2, byrow=T)

    c_ncross<-convex_hull_check(rpca, hline, t(cross_side))

    if(length(which(c_ncross==T)) %% 2 != 0){return(T)}
    else{return(F)}

  })

  return(exist)

}


#'Delete the same side in side set
#'@param sides a side set
#'@return a side set without overlapb
#'
#'
sides_set<-function(sides){

  if(!is.matrix(sides)){sides<-t(as.matrix(sides))}

  check_sides<-matrix(0, 2, ncol(sides)*2)

  t<-1

  for (i in 1:ncol(sides)) {
    check_sides[,t]<-c(sides[1,i], sides[2,i])
    t<-t+1
    check_sides[,t]<-c(sides[3,i], sides[4,i])
    t<-t+1
  }

  for(j in seq(ncol(check_sides), 1)){

    for(k in seq((ncol(check_sides)-1), 1)){

      if(j!=k && setequal(check_sides[,j], check_sides[,k])){

        check_sides<-check_sides[,-j]

        break

      }

    }

  }

  return(check_sides)

}


#'Distinguishing whether interpolated points are in a convex made of projected neighbor points via PCA from origin data.
#'@param rpca projected neighbor points via PCA
#'@param hline horizontal lines from interpolated points
#'@param sides a side set of a convex
#'@return True values indicate whether horizontal lines from interpolated points cross sides of convex in the side set
#'
#'
convex_hull_check<-function(rpca, hline, sides){

  t1<-sapply(sides[,1], function(side){

    return((hline[1,1]-hline[2,1])*(rpca[["x"]][side, 2]-hline[1,2]))

  })

  t2<-sapply(sides[,2], function(side){

    return((hline[1,1]-hline[2,1])*(rpca[["x"]][side, 2]-hline[1,2]))

  })


  t3<-sapply(1:nrow(sides), function(k){

    return((rpca[["x"]][sides[k,1], 1]-rpca[["x"]][sides[k,2], 1])*(hline[1,2]-rpca[["x"]][sides[k,1], 2])+(rpca[["x"]][sides[k,1], 2]-rpca[["x"]][sides[k,2], 2])*(rpca[["x"]][sides[k,1], 1]-hline[1,1]))

  })


  t4<-sapply(1:nrow(sides), function(k){

    return((rpca[["x"]][sides[k,1], 1]-rpca[["x"]][sides[k,2], 1])*(hline[2,2]-rpca[["x"]][sides[k,1], 2])+(rpca[["x"]][sides[k,1], 2]-rpca[["x"]][sides[k,2], 2])*(rpca[["x"]][sides[k,1], 1]-hline[2,1]))

  })

  ncross<-((t1*t2)<0 & (t3*t4)<0)

  return(ncross)

}

#'Projecting intepolated points to coordinate system of original data
#'@param rpca points whose dimension are reduced via PCA
#'@param incord interpolated points
#'@param ori_cen a center point of neighbors
#'@return coordinates of intepolated points in coordinate system of original data
#'
origin_coordinate<-function(rpca, incord, ori_cen){

  eigen01<-as.matrix(rpca$rotation[,1])
  eigen02<-as.matrix(rpca$rotation[,2])

  if(!is.matrix(incord)){incord<-t(as.matrix(incord))}

  oricord<-sapply(1:nrow(incord), function(l){

    return((incord[l, 1]*(eigen01)+incord[l, 2]*(eigen02))+ori_cen)

  })

  return(t(oricord))

}

#' Deriving sides having a specified vertex of a convex
#' @param chul vertex set of a convex
#' @param vertx a vertex to derive what sides have this vertex
#' @importFrom magrittr %>%
#' @return side set having "vertx"
#'
convex_hull_vertx<-function(chul, vertx){


  vert_set<-c(chul[length(chul)], chul, chul[1])

  #debugText(vert.set, which(chul==vertx))

  sides<-matrix(0, 2, 2)
  sides[1,]<-which(chul==vertx) %>% '+'(., c(0, 1)) %>%   vert_set[.]
  sides[2,]<-which(chul==vertx) %>% '+'(., c(1, 2)) %>%   vert_set[.]

  return(sides)

}
