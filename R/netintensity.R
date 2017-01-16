#' Intensity Calculation over networks
#' 
#' Calculates edgewise and mean nodewise intensity fucntion for undirected networks
#' @author Matthias Eckardt
#' @param adj A adjacency matrix encoding network structure
#' @param x.coord.node X-coordinates of node in network
#' @param y.coord.node Y-coordinates of node in network
#' @param x.event X-coordinates of event
#' @param y.event Y-coordinates of event
#' 
#' @return counts vector containing the mean.intensity per node
#' @return edge.counts array containing edge.intensity for adjacent edges
#' @return g underlying graph object  
#' @return deg Vector of degrees per node


netintensity <- function(adj, x.coord.node, y.coord.node, x.event, y.event){
  require(spatstat)
  require(igraph)
  require(intervals)
  a.1 <- min(as.numeric(x.coord.node))
  a.2 <- max(as.numeric(x.coord.node))
  b.1 <- min(as.numeric(y.coord.node))
  b.2 <- max(as.numeric(y.coord.node))
  g <- graph_from_adjacency_matrix(as.matrix(adj), mode="undirected")
  coord.x <- cbind(V(g),x.coord.node)
  rownames(coord.x)<-seq(1,length(V(g)),1)
  coord.y <- cbind(V(g),y.coord.node)
  events<-rep(1,length(x.event)) 
  dats<-cbind(x.event, y.event,events)
  event.loc.interval<-Intervals(matrix(cbind(x.event, y.event), ncol=2))
  interval.area <- pairdist(ppp(x.coord.node,y.coord.node,xrange=c(a.1,a.2),yrange=c(b.1,b.2)))
  ne.g<-list()
  edge.counts <- list()
  nach <- list()
  counts <- list()
  indicator <- list()
  res <- list()
  mean.counts <- list()
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
      mat.x <- mat.y <-matrix(0, nrow=length(c.x), ncol=2)
      # defining intervals depending on position 
      for(j in 1:length(c.x)){
        xval<-c.x[[j]]
        if(xval<node.x){
          mat.x[j,]<-cbind(xval, node.x)
        }else{mat.x[j,]<-cbind(node.x,xval)
        }
        yval<-c.y[[j]]
        if(yval<node.y){
          mat.y[j,]<-cbind(yval, node.y)
        }else{mat.y[j,]<-cbind(node.y,yval)
        }
      }
      interval.node.x <- Intervals(mat.x)
      interval.node.y <- Intervals(mat.y)
      for(k in 1:length(c.x))
      {
        nach.names <- as.numeric(id.nach[k])
        res[k]<-as.numeric(abs(interval.area[as.numeric(id.node),as.numeric(nach.names)]))
        res<-as.numeric(res[k])
        x.a <- mat.x[,1]
        x.b <- mat.x[,2]
        y.a <- mat.y[,1]
        y.b <- mat.y[,2]
        x.aa<-matrix(rep(x.a, each=length(events)), ncol=length(x.a))    # repeat edge limits as edges times dim(events)
        x.bb<-matrix(rep(x.b, each=length(events)), ncol=length(x.a))
        y.aa<-matrix(rep(y.a, each=length(events)), ncol=length(x.a))
        y.bb<-matrix(rep(y.b, each=length(events)), ncol=length(x.a))
        x1 <- x.aa[,k]
        x2 <- x.bb[,k]
        y1 <- y.aa[,k]
        y2 <- y.bb[,k]
        dat.x<-dats[,1]
        dat.y<-dats[,2]
        dats.cbind <- cbind(x1,x2,y1,y2,dat.x,dat.y,events)       
        ev.mat <- matrix(0,ncol=length(c.x))   
        colnames(ev.mat)<-id.nach                 
        for(l in 1:dim(dats)[1])
        {
          dx.l <- dats.cbind[,5]
          dy.l <- dats.cbind[,6]
          x.1 <-  dats.cbind[,1]
          x.2 <-  dats.cbind[,2]
          y.1 <-  dats.cbind[,3]
          y.2 <-  dats.cbind[,4]
          ev <-   dats.cbind[,7]
          indicator[[l]]<- ev[x.1 <= dx.l & x.2 >= dx.l & y.1 <= dy.l & y.2 >= dy.l]  
          if(length(indicator[[l]])==0){
            indicator[[l]]<- 0
          }                                                               
        }
        ev.mat[k]<-Reduce('+', indicator[[l]])/res 
      }
      edge.counts[[i]] <- ev.mat   
      counts[[i]] <- Reduce('+', ev.mat)/deg[i]
      if(length(counts[[i]])==0){
        counts[[i]]<-0
      }
    }
    if(isn==1){      # Counts for isolated vertices = 0
      counts[[i]]<-0
    }
    cat("Calculation for node", i, "\n")
  } 
  
  list(mean.intensity = counts, edge.intensity=edge.counts, graph=g, n.neighbors= deg)
}
