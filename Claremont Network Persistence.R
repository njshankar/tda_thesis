# ----------------------- 0: IMPORT LIBRARIES ----------------------------------

library(tidyverse)
library(TDA)
library(scatterplot3d)
library(scales)
library(gganimate)
library(gridExtra)
library(patchwork)
library(igraph)
library(magrittr)


setwd("/Users/njshankar/Documents/School/Senior Thesis/Computations")

# ----------------------- 1: READ-IN DATA --------------------------------------

coauth <- read_csv("coauthor.csv", col_names = FALSE)
cdist <- read_csv("colabdist.csv", col_names = FALSE)
labels <- read_csv("labels.csv")

# ----------------------- 2: TREAT DATA ----------------------------------------

# ----------------------- 2.1: Matrix Forms ------------------------------------

get_distmat <- function(){
  
  coauth_mat <- matrix(0, 47, 47)
  cdist_mat <- matrix(0, 47, 47)
  dist_mat <- matrix(0, 47, 47)
  
  for (i in 1:46){
    
    for (j in i:47){
      
      coauth_mat[i,j] <- as.numeric(coauth[i,j])
      cdist_mat[i,j] <- as.numeric(cdist[i,j])
      
      coauth_mat[j,i] <- coauth_mat[i,j]
      cdist_mat[j,i] <- cdist_mat[i,j]
      
      dist_mat[i,j] <- min(1/coauth_mat[i,j], cdist_mat[i,j])
      dist_mat[j,i] <- dist_mat[i,j]
      
    }
    
  }
  
  return(dist_mat)
  
}


# ----------------------- 2.2: Compute Homology --------------------------------

compute_hom <- function(mdim, mscale){
  
  data <- get_distmat()
  
  pd <- ripsDiag(X = data, maxdimension = mdim, 
                 maxscale = mscale, dist = "arbitrary", 
                 library = "Dionysus", location = TRUE)
  
  return(pd)
  
}


# ----------------------- 2.3: Plot Barcode ------------------------------------

plot_barcode <- function(mdim, mscale){
  
  pd <- compute_hom(mdim, mscale)
  
  plot(pd[["diagram"]], barcode = TRUE, 
       main = "Barcode for CCM Network", cex.main = 2)
  
}

# ----------------------- 2.4: Plot Persistence Diagram ------------------------

plot_pd <- function(mdim, mscale){
  
  pd <- compute_hom(mdim, mscale)
  
  plot(pd[["diagram"]], 
       main = "Persistence for CCM Network", cex.main = 2)
  
}


# ----------------------- 3: MAIN CODE -----------------------------------------

hom <- compute_hom(3, 8)
barcode <- plot_barcode(3, 8)
pd <- plot_pd(3, 8)

# ----------------------- 3.1: Members of Clusters -----------------------------

get_sig_clusters <- function(thresh){
  
  sig_clusters <- c()
  
  for(i in 1: nrow(hom[["diagram"]])){
    
    if((hom[["diagram"]][i,1] == 0) & (hom[["diagram"]][i,3] - 
                                       hom[["diagram"]][i,2] >= thresh)){
      
      sig_clusters <- c(sig_clusters, i)
      
    }
    
  }
  
  return(sig_clusters)

}


get_cluster_members <- function(thresh){
  
  sig_clusters <- get_sig_clusters(thresh)
  
  print(sig_clusters)
  
  ids <- c()
  members <- c()
  
  for (i in sig_clusters[2:length(sig_clusters)]){
    
    ids <- c(ids, rep(i, nrow(hom[["cycleLocation"]][[i]])))
    members <- c(members, hom[["cycleLocation"]][[i]][,1])
    
  }
  
  return(data.frame(ids = ids, members = labels[members, 2]))
  
}



# ----------------------- 4: PLOTTING AS GRAPH ---------------------------------

network_filt <- function(r, distmat){
  
  adjmat <- apply(distmat, c(1,2), function(n) as.numeric(n <= r))
  
  filt <- graph_from_adjacency_matrix(adjmatrix = adjmat, 
                                      mode = "undirected",
                                      weighted = NULL,
                                      diag = FALSE,
                                      add.colnames = NULL)
  
  return(filt)
  
}




D = get_distmat()

g0 <- network_filt(0, D)
g1 <- network_filt(0.5, D)
g2 <- network_filt(1, D)
g3 <- network_filt(2, D)
g4 <- network_filt(3, D)
g5 <- network_filt(4, D)
g6 <- network_filt(5, D)
g7 <- network_filt(6, D)
g8 <- network_filt(7, D)

l = layout_nicely(g5)



gplot <- function(g, r){
  
  plot.igraph(g, layout = layout_nicely(g), edge.width = 3, 
              edge.color = "cornflowerblue", vertex.color = "darkorange",
              vertex.frame.color = "darkorange",
              main = paste("colab network: r =", r),
              vertex.label = labels$brief, 
              vertex.size = 0, vertex.label.cex = 0.63, 
              vertex.label.dist = 0.50,
              vertex.label.color = "firebrick")
  
}


gplot(g0, 0)
gplot(g1, 0.5)
gplot(g2, 1)
gplot(g3, 2)
gplot(g4, 3)
gplot(g5, 4)
gplot(g6, 5)
gplot(g7, 6)
gplot(g8, 7)


