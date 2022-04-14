# ----------------------- 0: IMPORT LIBRARIES ----------------------------------

library(tidyverse)
library(TDA)
library(scatterplot3d)
library(scales)
library(gganimate)
library(gridExtra)
library(patchwork)


setwd("/Users/njshankar/Documents/School/Senior Thesis/Computations")

# ----------------------- 1: FUNCTIONS TO GENERATE DATA ------------------------

# These functions generate n random observations from the specified shapes, and
# return a dataframe of data points. 

# ----------------------- 1.1: Data from 2-Disk --------------------------------

draw_disk2 <- function(n, cx, cy){
  
  theta = runif(n, min = 0, max = 2*pi)
  rad = runif(n, min = 0, max = 1)
  x = rad*cos(theta) + cx
  y = rad*sin(theta) + cy
  
  return(data.frame(x = x, y = y))
  
}

# ----------------------- 1.2: Data from 1-Sphere ------------------------------

draw_sphere1 <- function(n, cx, cy, err){
  
  theta = runif(n, min = 0, max = 2*pi)
  rad = rnorm(n, mean = 1, sd = err)
  x = rad*cos(theta) + cx
  y = rad*sin(theta) + cy
  
  return(data.frame(x = x, y = y))
  
}

# ----------------------- 1.3: Data from 2-Sphere ------------------------------

draw_sphere2 <- function(n, cx, cy, cz, err){
  
  theta = runif(n, min = 0, max = 2*pi)
  phi = runif(n, min = -pi, max = pi)
  rho = rnorm(n, mean = 1, sd = err)
  x = rho*sin(phi)*cos(theta) + cx
  y = rho*sin(phi)*sin(theta) + cy
  z = rho*cos(phi) + cz
  
  return(data.frame(x = x, y = y, z = z))
  
}

# ----------------------- 1.4: Data from 2-Torus -------------------------------

draw_torus2 <- function(n, cx, cy, cz, err){
  
  theta = runif(n, min = 0, max = 2*pi)
  phi = runif(n, min = 0, max = 2*pi)
  R = rnorm(n, mean = 1, sd = err)
  r = rnorm(n, mean = 0.25, sd = err/4)
  
  x = (R + r*cos(theta))*cos(phi)
  y = (R + r*cos(theta))*sin(phi)
  z = sin(theta)
  
  return(data.frame(x = x, y = y, z = z))
  
}

# ----------------------- 2: FUNCTIONS TO GENERATE 2D PLOTS --------------------

# ----------------------- 2.1: Generate Baseplot -------------------------------

gen_baseplot <- function(point_df){
  
  p <- ggplot(data = point_df) + 
    geom_point(mapping = aes(x = x, y = y, alpha = 0.5)) + 
    coord_fixed() +
    theme_minimal() + 
    xlab("X") + 
    ylab("Y") + 
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5,
                                                              size = 2,
                                                              face = "bold")) + 
    xlim(-2.5, 2.5) + 
    ylim(-1.2, 1.2)
  
  return(p)
  
}

# ----------------------- 2.2: Generate Circular Outline -----------------------

gen_circoutline <- function(cx, cy){
  
  x <- cos(seq(from = 0, to = 2*pi, length.out = 100)) + cx
  y <- sin(seq(from = 0, to = 2*pi, length.out = 100)) + cy
  
  circoutline <- data.frame(x = x, y = y)
  
  return(circoutline)
  
}

# ----------------------- 3: FUNCTIONS TO GET PERSISTENCE ----------------------

# ----------------------- 3.1: Persistence on VR -------------------------------


get_persistence <- function(point_df, maxnum){
  
  pd <- ripsDiag(X = point_df, maxdimension = 2, maxscale = maxnum, 
             dist = "euclidean", library = c("GUDHI", "Dionysus"), 
             location = TRUE, printProgress = FALSE)
  
  return(pd)
  
}

# ----------------------- 3.2: Plot Barcode ------------------------------------

plot_barcode <- function(point_df, label, maxnum){
  
  pd <- get_persistence(point_df, maxnum)
  
  plot(pd[["diagram"]], barcode = TRUE, 
       main = paste("Barcode:", label), cex.main = 2)
  
}

# ----------------------- 3.3: Plot Persistence Diagram ------------------------

plot_pd <- function(point_df, label, maxnum){
  
  pd <- get_persistence(point_df, maxnum)
  
  plot(pd[["diagram"]], main = paste("Persistence:", label), cex.main = 2)

}

# ----------------------- 4: FUNCTIONS TO PLOT SKELETONS -----------------------

# ----------------------- 4.1: Compute Distance Matrix -------------------------

gen_distmatrix <- function(point_df){
  
  n <- nrow(point_df)
  
  distmat <- matrix(0,n,n)
  
  for (i in 1:(n-1)){
    
    for (j in (i+1):n){
      
      distmat[i,j] <- sqrt(
        (point_df[i,1] - point_df[j,1])^2 + (point_df[i,2] - point_df[j,2])^2)
      
      distmat[j,i] <- distmat[i,j]
      
    }
    
  }
  
  return(distmat)
  
}

# ----------------------- 4.2: Generate 1-Skeleton -----------------------------

get_1skeleton <- function(t, point_df){
  
  n = nrow(point_df)
  
  distmat = gen_distmatrix(point_df)

  x <- c()
  y <- c()
  xend <- c()
  yend <- c()
  
  for (i in 1:(n-1)){
    
    for (j in (i+1):n){
      
      if(distmat[i,j] <= t){
        
        x <- c(x, point_df[i,1])
        y <- c(y, point_df[i,2])
        xend <- c(xend, point_df[j,1])
        yend <- c(yend, point_df[j,2])
        
      }
      
    }
    
  }
  
  segpoints <- data.frame(x = x, y = y, xend = xend, yend = yend)
  
  return(segpoints)
  
}

# ----------------------- 4.2: Generate 2-Skeleton -----------------------------

get_2skeleton <- function(t, point_df){
  
  n = nrow(point_df)
  
  distmat = gen_distmatrix(point_df)
  
  id_tally <- 0
  id <- c()
  x <- c()
  y <- c()

  
  for (i in 1:(n-2)){
    
    for (j in (i+1):(n-1)){
      
      for (k in (j+1):n){
        
        if((distmat[i,j] <= t) & (distmat[j,k] <= t) & (distmat[i,k] <= t)){
          
          id_tally <- id_tally + 1
          id <- c(id, rep(id_tally, 3))
          x <- c(x, point_df[i,1], point_df[j,1], point_df[k,1])
          y <- c(y, point_df[i,2], point_df[j,2], point_df[k,2])
          
        }
      
      }
      
    }
    
  }
  
  
  polypoints <- data.frame(id = factor(id), x = x, y  = y)
  
  return(polypoints)
  
}


# ----------------------- 4.3: Plot 2-Skeleton ---------------------------------

plot_2skeleton <- function(t, point_df){
  
  oneskel <- get_1skeleton(t, point_df)
  twoskel <- get_2skeleton(t, point_df)
  
  
  p <- ggplot() + 
    geom_polygon(data = twoskel, 
                 mapping = aes(x = x, y = y, fill = "firebrick", group = id), 
                 size = 0) +
    geom_segment(data = oneskel, mapping = aes(x = x, y = y, 
                                                 xend = xend, yend = yend), 
                 colour = "seagreen4") + 
    geom_point(data = point_df, mapping = aes(x = x, y = y, alpha = 0.5)) + 
    coord_fixed() +
    theme_minimal() + 
    xlab("X") + 
    ylab("Y") + 
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5, 
                                                              size = 42,
                                                              face = "bold"))
  
  return(p)
  
}






# ----------------------- 5: MAIN CODE -----------------------------------------

n = 100

# ----------------------- 5.1: 2-Disk ------------------------------------------

disk2_data <- draw_disk2(n = n, cx = 0, cy = 0)

disk2_plot <- gen_baseplot(disk2_data) + 
  geom_path(data = gen_circoutline(cx = 0, cy = 0), 
            mapping = aes(x = x, y = y), colour = "blue", size = 4, alpha = 0.2) + 
  geom_polygon(data = gen_circoutline(cx = 0, cy = 0), 
               mapping = aes(x = x, y = y, alpha = 0.2)) + 
  ggtitle("2-Disk")

disk2_barcode <- function(){plot_barcode(disk2_data, "2-Disk", 1.5)}
disk2_pd <- function(){plot_pd(disk2_data, "2-Disk", 1.5)}

disk2_barcode()
disk2_pd()


disk2_twoskel1 <- plot_2skeleton(0.06, disk2_data) + 
  geom_path(data = gen_circoutline(cx = 0, cy = 0), 
            mapping = aes(x = x, y = y), colour = "blue", 
            size = 4, alpha = 0.2) + 
  geom_polygon(data = gen_circoutline(cx = 0, cy = 0), 
               mapping = aes(x = x, y = y, alpha = 0.2)) + 
  geom_point(data = disk2_data, mapping = aes(x = x, y = y, alpha = 0.5)) + 
  ggtitle("r = 0.06")

disk2_twoskel2 <- plot_2skeleton(0.12, disk2_data) + 
  geom_path(data = gen_circoutline(cx = 0, cy = 0), 
            mapping = aes(x = x, y = y), colour = "blue", 
            size = 4, alpha = 0.2) + 
  geom_polygon(data = gen_circoutline(cx = 0, cy = 0), 
               mapping = aes(x = x, y = y, alpha = 0.2)) + 
  geom_point(data = disk2_data, mapping = aes(x = x, y = y, alpha = 0.5)) + 
  ggtitle("r = 0.12")

disk2_twoskel3 <- plot_2skeleton(0.18, disk2_data) + 
  geom_path(data = gen_circoutline(cx = 0, cy = 0), 
            mapping = aes(x = x, y = y), colour = "blue", 
            size = 4, alpha = 0.2) +
  geom_polygon(data = gen_circoutline(cx = 0, cy = 0), 
               mapping = aes(x = x, y = y, alpha = 0.2)) + 
  geom_point(data = disk2_data, mapping = aes(x = x, y = y, alpha = 0.5)) + 
  ggtitle("r = 0.18")

disk2_twoskel4 <- plot_2skeleton(0.24, disk2_data) + 
  geom_path(data = gen_circoutline(cx = 0, cy = 0), 
            mapping = aes(x = x, y = y), colour = "blue", 
            size = 4, alpha = 0.2) +
  geom_polygon(data = gen_circoutline(cx = 0, cy = 0), 
               mapping = aes(x = x, y = y, alpha = 0.2)) + 
  geom_point(data = disk2_data, mapping = aes(x = x, y = y, alpha = 0.5)) + 
  ggtitle("r = 0.24")

disk2_twoskel5 <- plot_2skeleton(0.30, disk2_data) + 
  geom_path(data = gen_circoutline(cx = 0, cy = 0), 
            mapping = aes(x = x, y = y), colour = "blue", 
            size = 4, alpha = 0.2) +
  geom_polygon(data = gen_circoutline(cx = 0, cy = 0), 
               mapping = aes(x = x, y = y, alpha = 0.2)) + 
  geom_point(data = disk2_data, mapping = aes(x = x, y = y, alpha = 0.5)) + 
  ggtitle("r = 0.30")


disk2_twoskel6 <- plot_2skeleton(0.36, disk2_data) + 
  geom_path(data = gen_circoutline(cx = 0, cy = 0), 
            mapping = aes(x = x, y = y), colour = "blue", 
            size = 4, alpha = 0.2) +
  geom_polygon(data = gen_circoutline(cx = 0, cy = 0), 
               mapping = aes(x = x, y = y, alpha = 0.2)) + 
  geom_point(data = disk2_data, mapping = aes(x = x, y = y, alpha = 0.5)) + 
  ggtitle("r = 0.36")

disk2_twoskel7 <- plot_2skeleton(0.42, disk2_data) + 
  geom_path(data = gen_circoutline(cx = 0, cy = 0), 
            mapping = aes(x = x, y = y), colour = "blue", 
            size = 4, alpha = 0.2) +
  geom_polygon(data = gen_circoutline(cx = 0, cy = 0), 
               mapping = aes(x = x, y = y, alpha = 0.2)) + 
  geom_point(data = disk2_data, mapping = aes(x = x, y = y, alpha = 0.5)) + 
  ggtitle("r = 0.42")

disk2_twoskel8 <- plot_2skeleton(0.48, disk2_data) + 
  geom_path(data = gen_circoutline(cx = 0, cy = 0), 
            mapping = aes(x = x, y = y), colour = "blue", 
            size = 4, alpha = 0.2) +
  geom_polygon(data = gen_circoutline(cx = 0, cy = 0), 
               mapping = aes(x = x, y = y, alpha = 0.2)) + 
  geom_point(data = disk2_data, mapping = aes(x = x, y = y, alpha = 0.5)) + 
  ggtitle("r = 0.48")

disk2_plot
disk2_twoskel1
disk2_twoskel2
disk2_twoskel3
disk2_twoskel4
disk2_twoskel5
disk2_twoskel6
disk2_twoskel7
disk2_twoskel8


# ----------------------- 5.2: 1-Sphere ----------------------------------------

sphere1_data <- draw_sphere1(n = n, cx = 0, cy = 0, err = 0.05)

sphere1_plot <- gen_baseplot(sphere1_data) + 
  geom_path(data = gen_circoutline(cx = 0, cy = 0), 
            mapping = aes(x = x, y = y),  colour = "blue", 
            size = 4, alpha = 0.2) + 
  ggtitle("1-Sphere")

sphere1_barcode <- function(){plot_barcode(sphere1_data, "1-Sphere", 1.5)}
sphere1_pd <- function(){plot_pd(sphere1_data, "1-Sphere", 1.5)}

sphere1_barcode()
sphere1_pd()


sphere1_twoskel1 <- plot_2skeleton(0.06, sphere1_data) + 
  geom_point(data = sphere1_data, mapping = aes(x = x, y = y, alpha = 0.5)) +
  ggtitle("r = 0.06")

sphere1_twoskel2 <- plot_2skeleton(0.12, sphere1_data) + 
  geom_point(data = sphere1_data, mapping = aes(x = x, y = y, alpha = 0.5)) +
  ggtitle("r = 0.12")

sphere1_twoskel3 <- plot_2skeleton(0.18, sphere1_data) + 
  geom_point(data = sphere1_data, mapping = aes(x = x, y = y, alpha = 0.5)) +
  ggtitle("r = 0.18")

sphere1_twoskel4 <- plot_2skeleton(0.24, sphere1_data) + 
  geom_point(data = sphere1_data, mapping = aes(x = x, y = y, alpha = 0.5)) +
  ggtitle("r = 0.24")

sphere1_twoskel5 <- plot_2skeleton(0.30, sphere1_data) + 
  geom_point(data = sphere1_data, mapping = aes(x = x, y = y, alpha = 0.5)) +
  ggtitle("r = 0.30")


sphere1_twoskel6 <- plot_2skeleton(0.36, sphere1_data) + 
  geom_point(data = sphere1_data, mapping = aes(x = x, y = y, alpha = 0.5)) +
  ggtitle("r = 0.36")

sphere1_twoskel7 <- plot_2skeleton(0.42, sphere1_data) + 
  geom_point(data = sphere1_data, mapping = aes(x = x, y = y, alpha = 0.5)) +
  ggtitle("r = 0.42")

sphere1_twoskel8 <- plot_2skeleton(0.48, sphere1_data) + 
  geom_point(data = sphere1_data, mapping = aes(x = x, y = y, alpha = 0.5)) + 
  ggtitle("r = 0.48")

sphere1_plot
sphere1_twoskel1
sphere1_twoskel2
sphere1_twoskel3
sphere1_twoskel4
sphere1_twoskel5
sphere1_twoskel6
sphere1_twoskel7
sphere1_twoskel8

# ----------------------- 5.3: 2-Sphere ----------------------------------------

sphere2_data <- draw_sphere2(n = n, cx = 0, cy = 0, cz = 0, err = 0.05)

sphere2_barcode <- function(){plot_barcode(sphere2_data, "2-Sphere", 2)}
sphere2_pd <- function(){plot_pd(sphere2_data, "2-Sphere", 2)}

sphere2_barcode()
sphere2_pd()

# ----------------------- 5.4: 2-Torus -----------------------------------------

torus2_data <- draw_torus2(n = n, cx = 0, cy = 0, cz = 0, err = 0.05)

torus2_barcode <- function(){plot_barcode(torus2_data, "2-Torus", 2)}
torus2_pd <- function(){plot_pd(torus2_data, "2-Torus", 2)}

torus2_barcode()
torus2_pd()

# ----------------------- 5.5: Tangent Circles ---------------------------------

twocirc_data <- rbind(draw_sphere1(n = ceiling(n/2), cx = 0.5, 
                             cy = 0, err = 0.05),
                      draw_sphere1(n = ceiling(n/2), cx = -0.5, 
                                   cy = 0, err = 0.05))

twocirc_plot <- gen_baseplot(twocirc_data) + 
  geom_path(data = gen_circoutline(cx = 0.5, cy = 0), 
            mapping = aes(x = x, y = y), colour = "blue", 
            size = 4, alpha = 0.2) + 
  geom_path(data = gen_circoutline(cx = -0.5, cy = 0), 
            mapping = aes(x = x, y = y), colour = "blue", 
            size = 4, alpha = 0.2) + 
  ggtitle("Two Circles")


twocirc_barcode <- function(){plot_barcode(twocirc_data, "Two Circles", 1.5)}
twocirc_pd <- function(){plot_pd(twocirc_data, "Two Circles", 1.5)}

twocirc_barcode()
twocirc_pd()


twocirc_twoskel1 <- plot_2skeleton(0.12, twocirc_data) + 
  geom_point(data = twocirc_data, mapping = aes(x = x, y = y, alpha = 0.5)) +
  ggtitle("r = 0.12")

twocirc_twoskel2 <- plot_2skeleton(0.24, twocirc_data) + 
  geom_point(data = twocirc_data, mapping = aes(x = x, y = y, alpha = 0.5)) +
  ggtitle("r = 0.24")

twocirc_twoskel3 <- plot_2skeleton(0.36, twocirc_data) + 
  geom_point(data = twocirc_data, mapping = aes(x = x, y = y, alpha = 0.5)) +
  ggtitle("r = 0.36")

twocirc_twoskel4 <- plot_2skeleton(0.48, twocirc_data) + 
  geom_point(data = twocirc_data, mapping = aes(x = x, y = y, alpha = 0.5)) +
  ggtitle("r = 0.48")

twocirc_twoskel5 <- plot_2skeleton(0.60, twocirc_data) + 
  geom_point(data = twocirc_data, mapping = aes(x = x, y = y, alpha = 0.5)) +
  ggtitle("r = 0.60")


twocirc_twoskel6 <- plot_2skeleton(0.72, twocirc_data) + 
  geom_point(data = twocirc_data, mapping = aes(x = x, y = y, alpha = 0.5)) +
  ggtitle("r = 0.72")

twocirc_twoskel7 <- plot_2skeleton(0.84, twocirc_data) + 
  geom_point(data = twocirc_data, mapping = aes(x = x, y = y, alpha = 0.5)) +
  ggtitle("r = 0.84")

twocirc_twoskel8 <- plot_2skeleton(0.96, twocirc_data) + 
  geom_point(data = twocirc_data, mapping = aes(x = x, y = y, alpha = 0.5)) + 
  ggtitle("r = 0.96")

twocirc_plot
twocirc_twoskel1
twocirc_twoskel2
twocirc_twoskel3
twocirc_twoskel4
twocirc_twoskel5
twocirc_twoskel6
twocirc_twoskel7
twocirc_twoskel8

# ----------------------- 5.6: Separated Circles -------------------------------


sepcirc_data <- rbind(draw_sphere1(n = ceiling(n/2), cx = 1.5, 
                                   cy = 0, err = 0.05),
                      draw_sphere1(n = ceiling(n/2), cx = -1.5, 
                                   cy = 0, err = 0.05))

sepcirc_plot <- gen_baseplot(sepcirc_data) + 
  geom_path(data = gen_circoutline(cx = 1.5, cy = 0), 
            mapping = aes(x = x, y = y), colour = "blue", 
            size = 4, alpha = 0.2) + 
  geom_path(data = gen_circoutline(cx = -1.5, cy = 0), 
            mapping = aes(x = x, y = y), colour = "blue", 
            size = 4, alpha = 0.2) + 
  ggtitle("Separated Circles")


sepcirc_barcode <- function(){plot_barcode(sepcirc_data, "Separated Circles", 1.5)}
sepcirc_pd <- function(){plot_pd(sepcirc_data, "Separated Circles", 1.5)}

sepcirc_barcode()
sepcirc_pd()


sepcirc_twoskel1 <- plot_2skeleton(0.12, sepcirc_data) + 
  geom_point(data = sepcirc_data, mapping = aes(x = x, y = y, alpha = 0.5)) +
  ggtitle("r = 0.12")

sepcirc_twoskel2 <- plot_2skeleton(0.24, sepcirc_data) + 
  geom_point(data = sepcirc_data, mapping = aes(x = x, y = y, alpha = 0.5)) +
  ggtitle("r = 0.24")

sepcirc_twoskel3 <- plot_2skeleton(0.36, sepcirc_data) + 
  geom_point(data = sepcirc_data, mapping = aes(x = x, y = y, alpha = 0.5)) +
  ggtitle("r = 0.36")

sepcirc_twoskel4 <- plot_2skeleton(0.48, sepcirc_data) + 
  geom_point(data = sepcirc_data, mapping = aes(x = x, y = y, alpha = 0.5)) +
  ggtitle("r = 0.48")

sepcirc_twoskel5 <- plot_2skeleton(0.60, sepcirc_data) + 
  geom_point(data = sepcirc_data, mapping = aes(x = x, y = y, alpha = 0.5)) +
  ggtitle("r = 0.60")


sepcirc_twoskel6 <- plot_2skeleton(0.72, sepcirc_data) + 
  geom_point(data = sepcirc_data, mapping = aes(x = x, y = y, alpha = 0.5)) +
  ggtitle("r = 0.72")

sepcirc_twoskel7 <- plot_2skeleton(0.84, sepcirc_data) + 
  geom_point(data = sepcirc_data, mapping = aes(x = x, y = y, alpha = 0.5)) +
  ggtitle("r = 0.84")

sepcirc_twoskel8 <- plot_2skeleton(0.96, sepcirc_data) + 
  geom_point(data = sepcirc_data, mapping = aes(x = x, y = y, alpha = 0.5)) + 
  ggtitle("r = 0.96")

sepcirc_plot
sepcirc_twoskel1
sepcirc_twoskel2
sepcirc_twoskel3
sepcirc_twoskel4
sepcirc_twoskel5
sepcirc_twoskel6
sepcirc_twoskel7
sepcirc_twoskel8


# ----------------------- 5.7: Exporting Barcodes ------------------------------

tiff("barcode1.tiff", units="in", width=8, height=5, res=600)

plot_barcode(disk2_data, "2-Disk", 1.5)

dev.off()

tiff("barcode2.tiff", units="in", width=8, height=5, res=600)

plot_barcode(sphere1_data, "1-Sphere", 1.5)

dev.off()

tiff("barcode3.tiff", units="in", width=8, height=5, res=600)

plot_barcode(sphere2_data, "2-Sphere", 1.5)

dev.off()

tiff("barcode4.tiff", units="in", width=8, height=5, res=600)

plot_barcode(torus2_data, "2-Torus", 1.5)

dev.off()

tiff("barcode5.tiff", units="in", width=8, height=5, res=600)

plot_barcode(twocirc_data, "Two Circles", 1.5)

dev.off()

tiff("barcode6.tiff", units="in", width=8, height=5, res=600)

plot_barcode(sepcirc_data, "Separated Circles", 1.5)

dev.off()


# ----------------------- 5.8: Exporting PD's ----------------------------------

tiff("pd1.tiff", units="in", width=8, height=5, res=600)

plot_pd(disk2_data, "2-Disk", 1.5)

dev.off()

tiff("pd2.tiff", units="in", width=8, height=5, res=600)

plot_pd(sphere1_data, "1-Sphere", 1.5)

dev.off()

tiff("pd3.tiff", units="in", width=8, height=5, res=600)

plot_pd(sphere2_data, "2-Sphere", 1.5)

dev.off()

tiff("pd4.tiff", units="in", width=8, height=5, res=600)

plot_pd(torus2_data, "2-Torus", 1.5)

dev.off()

tiff("pd5.tiff", units="in", width=8, height=5, res=600)

plot_pd(twocirc_data, "Two Circles", 1.5)

dev.off()

tiff("pd6.tiff", units="in", width=8, height=5, res=600)

plot_pd(sepcirc_data, "Separated Circles", 1.5)

dev.off()


# ----------------------- 5.9: Exporting Sample Plots --------------------------

tiff("sp1.tiff", units="in", width=8, height=5, res=600)

disk2_plot

dev.off()

tiff("sp2.tiff", units="in", width=8, height=5, res=600)

sphere1_plot

dev.off()

tiff("sp3.tiff", units="in", width=8, height=5, res=600)

twocirc_plot

dev.off()

tiff("sp4.tiff", units="in", width=8, height=5, res=600)

sepcirc_plot

dev.off()

# ----------------------- 5.9: Exporting Skeletons -----------------------------

tiff("skel11.tiff", units="in", width=8, height=5, res=600)

disk2_twoskel1

dev.off()

tiff("skel12.tiff", units="in", width=8, height=5, res=600)

disk2_twoskel2

dev.off()

tiff("skel13.tiff", units="in", width=8, height=5, res=600)

disk2_twoskel3

dev.off()

tiff("skel14.tiff", units="in", width=8, height=5, res=600)

disk2_twoskel4

dev.off()

tiff("skel15.tiff", units="in", width=8, height=5, res=600)

disk2_twoskel5

dev.off()

tiff("skel16.tiff", units="in", width=8, height=5, res=600)

disk2_twoskel6

dev.off()

tiff("skel17.tiff", units="in", width=8, height=5, res=600)

disk2_twoskel7

dev.off()

tiff("skel18.tiff", units="in", width=8, height=5, res=600)

disk2_twoskel8

dev.off()







tiff("skel21.tiff", units="in", width=8, height=5, res=600)

sphere1_twoskel1

dev.off()

tiff("skel22.tiff", units="in", width=8, height=5, res=600)

sphere1_twoskel2

dev.off()

tiff("skel23.tiff", units="in", width=8, height=5, res=600)

sphere1_twoskel3

dev.off()

tiff("skel24.tiff", units="in", width=8, height=5, res=600)

sphere1_twoskel4

dev.off()

tiff("skel25.tiff", units="in", width=8, height=5, res=600)

sphere1_twoskel5

dev.off()

tiff("skel26.tiff", units="in", width=8, height=5, res=600)

sphere1_twoskel6

dev.off()

tiff("skel27.tiff", units="in", width=8, height=5, res=600)

sphere1_twoskel7

dev.off()

tiff("skel28.tiff", units="in", width=8, height=5, res=600)

sphere1_twoskel8

dev.off()








tiff("skel31.tiff", units="in", width=8, height=5, res=600)

twocirc_twoskel1

dev.off()

tiff("skel32.tiff", units="in", width=8, height=5, res=600)

twocirc_twoskel2

dev.off()

tiff("skel33.tiff", units="in", width=8, height=5, res=600)

twocirc_twoskel3

dev.off()

tiff("skel34.tiff", units="in", width=8, height=5, res=600)

twocirc_twoskel4

dev.off()

tiff("skel35.tiff", units="in", width=8, height=5, res=600)

twocirc_twoskel5

dev.off()

tiff("skel36.tiff", units="in", width=8, height=5, res=600)

twocirc_twoskel6

dev.off()

tiff("skel37.tiff", units="in", width=8, height=5, res=600)

twocirc_twoskel7

dev.off()

tiff("skel38.tiff", units="in", width=8, height=5, res=600)

twocirc_twoskel8

dev.off()








tiff("skel41.tiff", units="in", width=8, height=5, res=600)

sepcirc_twoskel1

dev.off()

tiff("skel42.tiff", units="in", width=8, height=5, res=600)

sepcirc_twoskel2

dev.off()

tiff("skel43.tiff", units="in", width=8, height=5, res=600)

sepcirc_twoskel3

dev.off()

tiff("skel44.tiff", units="in", width=8, height=5, res=600)

sepcirc_twoskel4

dev.off()

tiff("skel45.tiff", units="in", width=8, height=5, res=600)

sepcirc_twoskel5

dev.off()

tiff("skel46.tiff", units="in", width=8, height=5, res=600)

sepcirc_twoskel6

dev.off()

tiff("skel47.tiff", units="in", width=8, height=5, res=600)

sepcirc_twoskel7

dev.off()

tiff("skel48.tiff", units="in", width=8, height=5, res=600)

sepcirc_twoskel8

dev.off()