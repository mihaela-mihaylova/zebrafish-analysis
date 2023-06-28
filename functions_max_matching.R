library( celltrackR )
library(tidyverse)
library(igraph)

# Function that computes the distances between the points in all frames where they are visible together.
# It returns the minimum distance for the whole track (this is used later when we filter for > 10 micrometer distance)
t.dist <- function(t1, t2){ 
  x <- dist.trace( t1, t2 )
  if( length(x) == 0 ){
    return(Inf)
  }
  min(x)
}

# calculates Euclidean distance for the cases in which points from the two tracks are visible in the same frame
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
    return(0)
  } 
  
  # generate complete distance matrix between tracks
  # the distance between tracks is defined as Inf
  # if the tracks do not have points present in the same frame
  # and as the smallest distance between points from the two tracks
  # which are in the same frame
  xx <- sapply( t2, function(x) sapply(t1, function(y) t.dist( x, y ) ))
  ## Drop cells that are not closer than min_dist_cells pixels to any other cell in
  ## the other set
  if(is.null(dim(xx))|all(is.infinite(xx))){
    return(0)
  }
  t1 <- t1[apply( xx, 1, min ) < min_dist_cells]
  xx <- xx[apply( xx, 1, min ) < min_dist_cells,]

    if(is.null(dim(xx))|all(is.infinite(xx))){
    return(0)
  }
 
  t2 <- t2[apply( xx, 2, min ) < min_dist_cells]
  xx <- xx[,apply( xx, 2, min ) < min_dist_cells]
  c <- xx
  
  if(is.null(dim(xx))|all(is.infinite(xx))){
    return(0)
  }
  ## Focus on less than 10 pixels distances
  xx[xx>min_dist_cells] <- Inf
  if(is.null(dim(xx))|all(is.infinite(xx))){
    return(0)
  }
  # declare a matrix with the dimensions of xx but only with NA values
  mm <- matrix(NA,nrow(xx),ncol(xx))
  
  # extract only the positions (row, col) of entries in xx which have a finite (not Inf) value
  ii <- which( is.finite(xx), arr.ind=TRUE )
  if(is_empty(ii)){
    return(0)
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
  # generate a matrix with interaction lengths,
  # turn into a weighted bipartite graph with interaction lengths as weights
  # returns matching_weight/matching_size (the average edge weight in the maximum matching)
  
  
  # Subtract 1 from all interaction lengths:
  # - if we only have 1 frame, it is not clear if the cells interacted or just
  # happened to be close when the frame was shot
  # - for interactions of 2+ frames, we use the very first frame as an indication that the
  # interaction has started.
  
  # Multiply by 1.5 (min), as the video resolution is 90s
  mat <- (mm-1)*1.5
  mat[is.na(mm)] <- 0
 
  
  d <- dim(mat)
  numrow <-d[1]
  numcol <-d[2]
  row_names <- list()
  col_names <- list()
  
  for (i in 1:numrow) {
    row_names <- append(row_names, paste0('r',i))
  }
  
  for (i in 1:numcol) {
    col_names <- append(col_names, paste0('c',i))
  }
  
  colnames(mat) <- col_names
  rownames(mat) <- row_names
  
  tdf <- as.matrix(mat)
  # generate a bipartite graph from interaction duration matrix
  # the weights are the interaction durations
  g <- graph.incidence(tdf, weighted = TRUE, directed=FALSE)
  max_matching <- max_bipartite_match(g)
  # in case there are no cells close to each other
  if(max_matching$matching_weight==0){
    avg_edge_weight_in_max_matching <- 0
  } else{
    avg_edge_weight_in_max_matching <- max_matching$matching_weight/max_matching$matching_size
  }
  return(avg_edge_weight_in_max_matching)
}

# a function that takes in tibbles for the time intervals of interest, merges them in one df
# and returns the df
merge_interaction_data <- function(plt_02h, plt_24h, plt_46h, plt_68h){
  
  # check if a tibble is empty (meaning no two cells neut and macr are < 10micrometers apart)
  if(nrow(plt_02h)==0){
    #print(plt_02h)
    save_02h <- c(NA)
  } else{
    #print(plt_02h)
    save_02h <- list(plt_02h[[1]])[[1]]
  }
  
  if(nrow(plt_24h)==0){
    #print(plt_24h)
    save_24h <- c(NA)
  } else{
    #print(plt_24h)
    save_24h <- list(plt_24h[[1]])[[1]]
  }
  
  if(nrow(plt_46h)==0){
    #print(plt_46h)
    save_46h <- c(NA)
  } else{
    #print(plt_46h)
    save_46h <- list(plt_46h[[1]])[[1]]
  }
  
  if(nrow(plt_68h)==0){
    #print(plt_68h)
    save_68h <- c(NA)
  } else{
    #print(plt_68h)
    save_68h <- list(plt_68h[[1]])[[1]]
  }
  
  # find max length
  max_length <- max(length(save_02h), length(save_24h), length(save_46h), length(save_68h))
  
  # assign equal lengths
  length(save_02h) <- max_length
  length(save_24h) <- max_length
  length(save_46h) <- max_length
  length(save_68h) <- max_length
  all_interactions_per_int <- cbind(save_02h, save_24h, save_46h, save_68h )
  colnames(all_interactions_per_int) <- c("int_02h","int_24h", "int_46h","int_68h")
  return(all_interactions_per_int)
}