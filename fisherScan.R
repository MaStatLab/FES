fisherScan = function(x,y,K1=7,K2=7,fixed.M=FALSE,M=NULL,rank.transform=TRUE,plot.it=FALSE,plot.sidek=FALSE,plot.scale=FALSE,output.to.files=FALSE,mid.p=TRUE,compute.bf=FALSE,verbose=FALSE) {
  
  ## Preprocessing: marginally rank transform the data
 
  n=length(x)
  
  if (rank.transform) {
    ecdf.x = ecdf(x)
    ecdf.y = ecdf(y)
    x.tran = ecdf.x(x)-1/n
    y.tran = ecdf.y(y)-1/n
  }
  else {
    range.x = range(x)
    range.y = range(y)
    x.tran = (x-range.x[1])/(range.x[2]-range.x[1])*(1-1E-7)
    y.tran = (y-range.y[1])/(range.y[2]-range.y[1])*(1-1E-7)
  }
  
  ## Transform the data into a 2^K1 by 2^K2 contingency table
  n.list = list()
  for (k1 in 0:K1) {
    for (k2 in 0:K2) {
      index = make.index(k1,k2,K2)
      n.mat = matrix(0,2^k1,2^k2)
      
      for (obs in 1:n) {
        x.curr = x.tran[obs]
        y.curr = y.tran[obs]
        i = floor(x.curr*2^k1) + 1
        j = floor(y.curr*2^k2) + 1
        n.mat[i,j] = n.mat[i,j] + 1
      }
      
      n.list[[index]] = n.mat
    }
  }
  
  ## Start Fisher's exact scanning
  pval.list = list()
  if (compute.bf) bf.list = list()
  
  for (k1 in 0:(K1-1)) {
    for (k2 in 0:(K2-1)) {
      
      index.n = make.index(k1,k2,K2)
      n.mat.k1plus1.k2plus1 = n.list[[index.n + (K2+1) + 1]]
      pval.mat = matrix(NA,2^k1,2^k2)
      if (compute.bf) bf.mat = matrix(NA,2^k1,2^k2)
      
      for (i in 1:(2^k1)) {
        for (j in 1:(2^k2)) {
          twobytwotable = n.mat.k1plus1.k2plus1[c(2*i-1,2*i),c(2*j-1,2*j)]
          nx.tots = apply(twobytwotable,2,sum)
          ny.tots = apply(twobytwotable,1,sum)
          
          if (sum(twobytwotable)>25 && sum(nx.tots > 10)==2 && sum(ny.tots >10) == 2) {
            pval.mat[i,j] = fisher.test(twobytwotable,conf.int=FALSE)$p.value
            if (mid.p) pval.mat[i,j] = pval.mat[i,j] - dhyper(twobytwotable[1,1],nx.tots[1],nx.tots[2],ny.tots[1])/2
            
            if (compute.bf) bf.mat[i,j] = contingencyTableBF(twobytwotable,sampleType="hypergeom",priorConcentration=1)@bayesFactor$bf
          }
        }
      }
      
      index.p = make.index(k1,k2,K2-1)
      pval.list[[index.p]] = pval.mat
      names(pval.list)[index.p] = paste("(",k1,",",k2,")",sep="")
      if (compute.bf) {
        bf.list[[index.p]] = bf.mat
        names(bf.list)[index.p] = paste("(",k1,",",k2,")",sep="")
      }
    }
  }
  
  # first column is the minimum p-value in each stratum, second column is L(i,j) the windows tested in each stratum
  pvals.mat = cbind(unlist(lapply(pval.list,function(x){min(x,na.rm=TRUE)})),unlist(lapply(pval.list,function(x){sum(!is.na(x))})))
  
  # Three level Sidek p-values
  ## Per stratum
  p.sidek.stratum = 1 - (1-pvals.mat[,1])^pvals.mat[,2]
  p.sidek.stratum[pvals.mat[,1]==Inf] = Inf
  p.bonferroni.stratum = pvals.mat[,1]*pvals.mat[,2]
  p.bonferroni.stratum[pvals.mat[,1]==Inf] = Inf
  p.min.stratum = pvals.mat[,1]
  
  if (is.null(M)) R = K1+K2-1
  else R = min(K1+K2-1,M+1)
  
  T = double(R)
  
  
  ## Per resolution
  p.sidek.resol = double(R)
  p.bonferroni.resol = double(R)
  p.min.resol = double(R)
  for (r in 0:(R-1)) {
    p.sidek.resol.curr = Inf
    p.bonferroni.resol.curr = Inf
    p.min.resol.curr = Inf
    
    for (i in max(0,r-(K2-1)):min(K1-1,r)) {
      j = r - i
      stratum.curr = paste("(",i,",",j,")",sep="")
      p.sidek.resol.curr = min(p.sidek.resol.curr,p.sidek.stratum[stratum.curr])
      p.bonferroni.resol.curr = min(p.bonferroni.resol.curr,p.bonferroni.stratum[stratum.curr])
      p.min.resol.curr = min(p.min.resol.curr, p.min.stratum[stratum.curr])
      
      if (min(p.sidek.stratum[stratum.curr]) != Inf) {
        T[r+1] = T[r+1] + 1
      }
    }
    if (T[r+1]>0) {
      p.sidek.resol[r+1] = 1-(1-p.sidek.resol.curr)^(T[r+1])
      p.bonferroni.resol[r+1] = p.bonferroni.resol.curr*T[r+1]
      p.min.resol[r+1] = p.min.resol.curr
    }
    else {
      p.sidek.resol[r+1] = 1
      p.bonferroni.resol[r+1] = 1
      p.min.resol[r+1] = 1
    }
  }
  names(p.sidek.resol) = names(p.bonferroni.resol) = names(p.min.resol) = 0:(R-1)
  
  if (verbose) {
    print(cbind(p.sidek.resol,p.bonferroni.resol,p.min.resol,T))
  }
  ## Global
  if (!fixed.M) R1 = sum(T>0)
  else R1 = min(R,max(T>0))
  
  p.sidek.global = 1 - (1 - min(p.sidek.resol))^R1
  p.bonferroni.global = min(p.bonferroni.resol)*R1
  p.min.global = min(p.min.resol)
  
  
  ### The Meta p-value
  
  if (R1 > 0) {
    p.meta.global = pnorm(sum(qnorm(p.sidek.resol[T>0]))/sqrt(R1))
  } else {
    p.meta.global = 1
  }
  
  
  pvals.global = c(p.sidek.global,p.bonferroni.global,p.meta.global,p.min.global)
  names(pvals.global) = c("Sidek","Bonferroni","Meta","Min")
  if (verbose) print(pvals.global)
  
  
  if (plot.it) {
    if (plot.scale) layout(matrix(c(1,2),nrow=1,ncol=2),width=c(5,1),height=5)
    if (!compute.bf) test.list = pval.list
    else test.list = bf.list
    
    if (verbose) print(test.list)
    
    zlim=c(0,min(100,max(c(-log(0.0001)+1,-log(unlist(pval.list)),na.rm=TRUE),na.rm=TRUE)))
    if (plot.sidek) zlim=c(0,min(100,max(c(-log(0.0001)+1,-log(p.sidek.global),na.rm=TRUE),na.rm=TRUE)))
    # print(zlim)
    for (level in 0:(R1-1)) {  ##min((K1+K2-2),6)) {
      for (k1 in max(0,level - (K2-1)):min(level,K1-1)) {
        k2 = level - k1
        
        index.p = make.index(k1,k2,K2-1)
        
        if (sum(!is.na(test.list[[index.p]]))>0) {
          if (compute.bf) z = test.list[[index.p]]
          else z = -log(test.list[[index.p]])
          
          alpha.i.j = 1 - (1-0.05)^(1/R1*1/T[k1+k2+1]*1/sum(!is.na(test.list[[index.p]])))
          signif.ind = which(test.list[[index.p]] < alpha.i.j,arr.ind=TRUE)
          # z = pmin(z,100)
          pval.adj = 1 - (test.list[[index.p]])^(R1*T[k1+k2+1]*sum(!is.na(test.list[[index.p]])))
          
          if (plot.sidek) z= -log(pval.adj)
          
          if (output.to.files) {
            width=5
            height=5.5
            pdf(file=paste("Figures/k1=",k1,",k2=",k2,".pdf",sep=""),width=width,height=height)
            main=""
          } else {main=paste(names(test.list)[index.p],sep="")}
          col = topo.colors(100)
          par(mar=c(2,2,2,2))
          my.image(x=((1:2^k1)-0.5)/2^k1,y=((1:2^k2)-0.5)/2^k2,z=z,xlim=c(0,1),ylim=c(0,1),zlim=zlim,
                   xlab="X",ylab="Y",main=main,bf=compute.bf,plot.scale=plot.scale,col=col)
          grid(nx=2^k1,ny=2^k2,lty="solid")
          if (nrow(signif.ind)>=1) {
            par(new=TRUE)
            rect(xleft=(signif.ind[,1]-1)/2^k1,xright=signif.ind[,1]/2^k1,
                 ybottom=(signif.ind[,2]-1)/2^k2,ytop=signif.ind[,2]/2^k2,border="red",lwd=3)
          }
          if (output.to.files) dev.off()
          
        }
        
      }
    }
    if (output.to.files) pdf(file=paste("pval_scale.pdf",sep=""),width=height/2,height=height)
    
    par(mar=c(2,8,2,2))
    zlim=c(0,min(100,max(c(-log(0.0001)+1,-log(unlist(pval.list)),na.rm=TRUE),na.rm=TRUE)))
    if (plot.sidek) zlim=c(0,min(100,max(c(-log(0.0001)+1,-log(p.sidek.global),na.rm=TRUE),na.rm=TRUE)))
    
    image.scale(-log(unlist(pval.list)), zlim=zlim, col=col,axis.pos=2,add.axis=FALSE)
    
    ticks = 10^(seq(1,zlim[2]/log(10)))     #c(0.1,0.01,0.001,0.0001,0.00001,1e-6))
    axis(2,at=log(ticks),labels=1/ticks,las=2)
    
    if (output.to.files) dev.off()
  }
  
  return(pvals.global)
}





