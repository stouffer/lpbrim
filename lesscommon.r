
netToAdjacency = function(net,useSparse=FALSE){
    x <- as.factor(net[,1])
    y <- as.factor(net[,2])
    z <- net[,3]
    if(!useSparse){
    	 adj <- matrix(0,
    	               nrow=nlevels(x),
    	               ncol=nlevels(y),
    	               dimnames=list(levels(x),levels(y))
    	              )
    	 adj[cbind(x, y)] <- z
    }else{
    	adj <- sparseMatrix(i=as.integer(x),
    						j=as.integer(y),
    						x=z,
    						dims=c(nlevels(x),nlevels(y)),
                    	    dimnames=list(levels(x), levels(y))
                       	   )
    }
    adj
}

modulesToSmat = function(modules){
    comms <- unique(modules)
    nc <- length(comms)
    Smat <- matrix(0,nrow=length(names(modules)),ncol=nc,dimnames=list(names(modules),comms))
    Smat[cbind(names(modules),modules)] <- 1
    Smat
}

