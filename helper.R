make.index = function(k1,k2,m) {
  return(k1*(m+1) + k2 + 1)
}

proportion.leq.noNA = function(x,a=0.05) {
  mean(x <= a,na.rm=TRUE)
}

mymine=function(x,y,tempfile="competitors/test.csv",javafile="competitors/MINE.jar") {
  xx=cbind(x,y)
  write("x,y",file=tempfile)
  write(t(xx),sep=",",file=tempfile,ncol=2,append=T)
  command <- paste("java -jar",javafile,tempfile,"-allPairs",sep=' ')
  system(command)
  res=scan(paste(tempfile,",allpairs,cv=0.0,B=n^0.6,Results.csv",sep=""),what="",sep=",")
  val=as.numeric(res[11])
  return(val)
}

plot.roc = function(p, p0, fpr = seq(0,1,by=0.01), col="black", xlim=c(0,1), ylim=c(0,1), main="", lty=1,lwd=2) { ## helper function to plot ROC
  plot(fpr,ecdf(p)(quantile(p0,fpr,na.rm=TRUE)), type="l", xlab="False postive rate", ylab="True positive rate", col=col, xlim=xlim, ylim=ylim, main=main, lty=lty,lwd=lwd)
}

# The code for function "image.scale" including all comments is obtained from
# https://gist.github.com/menugget/7689145#file-image-scale-2-r

image.scale <- function(z, zlim, col = heat.colors(12),
                        breaks, axis.pos=1, add.axis=TRUE, ...){
  
  #This function creates a color scale for use with the image()
  #function. Input parameters should be consistent with those
  #used in the corresponding image plot. The "axis.pos" argument
  #defines the side of the axis. The "add.axis" argument defines
  #whether the axis is added (default: TRUE)or not (FALSE).
  
  
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  if(axis.pos %in% c(1,3)){ylim<-c(0,1); xlim<-range(breaks)}
  if(axis.pos %in% c(2,4)){ylim<-range(breaks); xlim<-c(0,1)}
  plot(1,1,t="n",ylim=ylim, xlim=xlim, axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)){
    if(axis.pos %in% c(1,3)){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(axis.pos %in% c(2,4)){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
  box()
  if(add.axis) {axis(axis.pos)}
}

my.image = function(x,y,z,col=topo.colors(100),zlim=NULL,bf=FALSE,plot.scale=TRUE,...) {
  if(plot.scale) {
    layout(matrix(c(1,2),nrow=1,ncol=2),width=c(5,1),height=5)
    par(mar=c(2,2,2,2))
  }
  
  if (is.null(zlim)) {
    if (bf) { # For Bayes factors
      zlim = range(z,na.rm=TRUE)
      zlim[1] = min(zlim[1],log(0.1))
      zlim[2] = max(zlim[2],log(10))
    }
    else { # For p-values
      zlim=c(0,min(100,max(-log(0.0001)+1,z,na.rm=TRUE)))
    }
  }
  
  image(x=x,y=y,z=z,zlim=zlim,col=col,...)
  if (plot.scale) {
    par(mar=c(2,2,2,2))
    image.scale(z, zlim=zlim, col=col,axis.pos=2,add.axis=FALSE)
    
    if (bf) { # For Bayes factors
      ticks = log(c(0.1,0.2,0.5,1,2,5,10))
      axis(2,at=ticks,labels=exp(ticks),las=2)
    }
    else { # For p-values
      ticks = 10^(seq(1,zlim[2]/log(10)))    
      axis(2,at=log(ticks),labels=1/ticks,las=2)
      
    }
  }
}



