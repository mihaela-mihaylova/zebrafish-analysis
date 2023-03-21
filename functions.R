library( celltrackR )
library(tidyverse)

# Function that computes the distances between the points in all frames where they are visible together.
t.dist <- function(t1, t2){ 
  x <- dist.trace( t1, t2 )
  if( length(x) == 0 ){
    return(Inf)
  }
  min(x)
}

dist.trace <- function(t1, t2){ 
  ttr <- merge(t1, t2, by="t")
  if( nrow(ttr) == 0 ){
    return( numeric(0) )
  } else {
    sqrt(rowSums((ttr[,c("x.x","y.x")]-ttr[,c("x.y","y.y")])^2))
  }
}

roundUp <- function(x,m) m*ceiling(x / m)

# A function which extracts only the tracks within certain time interval
# and then calculates only cell interactions within that interval.
# The "interaction" two cells is defined as the longest stretch of 
# consecutive frames in which the cells are at a distance of < 10 units.
extract_time_interval <- function(d1, d2, min_time, max_time, min_time_incl, max_time_incl, min_track_length, min_dist_cells){
  # extract only the detections which fall into the interval of interest
  if(min_time_incl){
    d1 <- d1[d1$POSITION_T >= min_time,]
    d2 <- d2[d2$POSITION_T >= min_time,]
  } else{
    d1 <- d1[d1$POSITION_T > min_time,]
    d2 <- d2[d2$POSITION_T > min_time,]
  }
  
  if(max_time_incl){
    d1 <- d1[d1$POSITION_T <= max_time,]
    d2 <- d2[d2$POSITION_T <= max_time,]
  } else{
    d1 <- d1[d1$POSITION_T < max_time,]
    d2 <- d2[d2$POSITION_T < max_time,]
  }
  
  
  t1 <- as.tracks( d1, id.column="TRACK_ID", time.column="POSITION_T",
                   pos.columns=c("POSITION_X","POSITION_Y") )	
  t2 <- as.tracks( d2, id.column="TRACK_ID", time.column="POSITION_T",
                   pos.columns=c("POSITION_X","POSITION_Y") )
  
  ## Remove short tracks
  t1 <- t1[sapply(t1,nrow)>min_track_length]
  #print(length(t1))
  t2 <- t2[sapply(t2,nrow)>min_track_length]
  #print(length(t2))

  if(length(t1)==0|length(t2)==0){
 #   print("stop")
    return(tibble(c()))
  } 
  
  # generate complete distance matrix between tracks
  # the distance between tracks is defined as Inf
  # if the tracks do not have points present in the same frame
  # and as the smallest distance between points from the two tracks
  # which are in the same frame
  xx <- sapply( t2, function(x) sapply(t1, function(y) t.dist( x, y ) ))
  #  View(xx)
  ## Drop cells that are not closer than min_dist_cells pixels to any other cell in
  ## the other set
  if(is.null(dim(xx))|all(is.infinite(xx))){
    return(tibble(c()))
  }
  t1 <- t1[apply( xx, 1, min ) < min_dist_cells]
  xx <- xx[apply( xx, 1, min ) < min_dist_cells,]
  if(is.null(dim(xx))|all(is.infinite(xx))){
    return(tibble(c()))
  }
  t2 <- t2[apply( xx, 2, min ) < min_dist_cells]
  xx <- xx[,apply( xx, 2, min ) < min_dist_cells]

  if(is.null(dim(xx))|all(is.infinite(xx))){
    return(tibble(c()))
  }
  ## Focus on less than 10 pixels distances
  xx[xx>min_dist_cells] <- Inf
  if(is.null(dim(xx))|all(is.infinite(xx))){
    return(tibble(c()))
  }
  # declare a matrix with the dimensions of xx but only with NA values
  mm <- matrix(NA,nrow(xx),ncol(xx))

  # extract only the positions (row, col) of entries in xx which have a finite (not Inf) value
  ii <- which( is.finite(xx), arr.ind=TRUE )
  if(is_empty(ii)){
    return(tibble(c()))
  }

  for( i in 1:nrow(ii) ){
    # find the actual distance trace between t1 and t2
    d <- dist.trace( t1[[ii[i,1]]], t2[[ii[i,2]]] )
    # assign NA to values greater than 10
    d[d>min_dist_cells] <- NA
    # find the longest consecutive stretch of non-NA values
    # which is the stretch in which the two tracks are < 10 pixels apart
    mm[ii[i,1], ii[i,2]] <- length(na.contiguous(d))
  }
  ii <- which(mm==max(mm[!is.na(mm)]), arr.ind=TRUE)
  interactions_count <- mm[mm>0]
  interactions_count <- tibble(int_length=unlist(interactions_count[!is.na(interactions_count)])) 
}