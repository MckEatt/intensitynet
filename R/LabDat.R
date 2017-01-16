#' Calculation of edge- and nodewise means for planar networks
#' 
#' Calculates edgewise and nodewise means of covariates for planar networks
#' @author Matthias Eckardt
#' @param adj A adjacency matrix encoding network structure
#' @param x.coord.node X-coordinates of node in network
#' @param y.coord.node Y-coordinates of node in network
#' @param x.event X-coordinates of event
#' @param y.event Y-coordinates of event
#' @param codats 
#' 
#' @return dats vector containing the edgewise means
#' @return datmeans vector containing the nodewise means per neighbours


LabDat <- function(adj, x.coord.node, y.coord.node, x.event, y.event, codats){
  require(spatstat)
  require(igraph)
  g <- graph_from_adjacency_matrix(as.matrix(adj), mode="undirected")
  coord.x <- cbind(V(g),x.coord.node)
  rownames(coord.x)<-seq(1,length(V(g)),1)
  coord.y <- cbind(V(g),y.coord.node)
  events<-rep(1,length(x.event)) # associate each coordinates as given events
  dats<-cbind(x.event, y.event,events, codats)
  # obtain set of neigbours for each vertex from adjmatrix
  ne.g <-list()
  id.ev <- list()
  dat.l <- list()
  id <- list()
  id.n <- list()
  edgesid <- list()
  edge.daten <- list()
  mat.x <- list()
  mat.y <- list()
  nach <- list()
  dat.mean <- matrix(ncol=ncol(codats), nrow=length(V(g)))
  dat.k <-list()
  x.aa <- list()
  x.bb <- list()
  y.aa <- list()
  y.bb <- list()
  gdat <- list()
  res <- list()
  dat.l <- list()
  deg <- degree(g)
  isolated.nodes <- as.list(ifelse(deg==0,1,0))
  for(i in 1:length(V(g)))
  {
    isn <- isolated.nodes[[i]]
    if(isn == 0){
      ne.g<-neighbors(g,i)
      id.nach<-as.vector(ne.g)
      id.node<-as.vector(i)
      lnach <- length(ne.g)
      c.xd <- matrix(coord.x[which(coord.x%in%id.nach),], ncol=2)
      c.yd <- matrix(coord.y[which(coord.y%in%id.nach),], ncol=2)
      if(dim(c.xd)[2]>1){
        c.x <- c.xd[,2]
        c.y <- c.yd[,2]
      }else{
        c.x <- c.xd[2]
        c.y <- c.yd[2]
      }
      node.xd <- coord.x[which(coord.x%in%id.node),]
      node.yd <- coord.y[which(coord.y%in%id.node),]
      node.x <- node.xd[[2]]
      node.y <- node.yd[[2]]
      mat.x[[i]] <- mat.y[[i]] <-matrix(0, nrow=length(c.x), ncol=2)
      # defining intervals depending on position
      for(j in 1:length(c.x)){
        xval<-c.x[[j]]
        if(xval<node.x){
          mat.x[[i]][j,]<-cbind(xval, node.x)
        }else{mat.x[[i]][j,]<-cbind(node.x,xval)
        }
        yval<-c.y[[j]]
        if(yval<node.y){
          mat.y[[i]][j,]<-cbind(yval, node.y)
        }else{mat.y[[i]][j,]<-cbind(node.y,yval)
        }
      }
      x.aa[[i]]<-matrix(rep(mat.x[[i]][,1], each=length(events)), ncol=length(mat.x[[i]][,1]))    # repeat edge limits as edges times dim(events)
      x.bb[[i]]<-matrix(rep(mat.x[[i]][,2], each=length(events)), ncol=length(mat.x[[i]][,1]))
      y.aa[[i]]<-matrix(rep(mat.y[[i]][,1], each=length(events)), ncol=length(mat.x[[i]][,1]))
      y.bb[[i]]<-matrix(rep(mat.y[[i]][,2], each=length(events)), ncol=length(mat.x[[i]][,1]))
      x1 <- as.vector(x.aa[[i]])
      x2 <- as.vector(x.bb[[i]])
      y1 <- as.vector(y.aa[[i]])
      y2 <- as.vector(y.bb[[i]])
      dat.x <- dats[,1]
      dat.y <- dats[,2]
      covs <- dats[,4:ncol(dats)]
      dat.k[[i]] <- cbind(x1,dat.x, x2,y1, dat.y, y2,events, covs)
      esta <- dat.k[[i]]
      gdat[[i]] <- esta[x1 <= dat.x & x2 >= dat.x & y1 <= dat.y & y2 >= dat.y, -c(1:7)]
    }
    if(isn==1){      # Counts for isolated vertices = 0
      gdat[[i]]<-matrix(0,ncol=ncol(codats), nrow=nrow(dats))
    }
    dat.mean[i,] <- colMeans(gdat[[i]])
    dat.mean[dat.mean=="NaN"] <- 0
    cat("Calculation for node", i, "\n")
  }
  list(dats=gdat, datmeans = dat.mean)
}
