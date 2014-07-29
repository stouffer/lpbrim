library(nnet)
library(snow)
library(snowfall)

source('common.r')
source('lesscommon.r')

bQW = function(adj,Smat)
{
   # calcluate the column sums in the interaction matrix for use later
   cs <- colSums(adj)
   
   # calculate the positive contribution to modularity of matching any two rows
   # also ignore the diagonal since it doesn't influence anything
   A <- (adj %*% t(adj)) / sum(cs*(cs-1))
   diag(A) <- 0
   
   # calculate the expected contribution to modularity that would occur if any two rows show overlapped by chance
   # also ignore the diagonal since it doesn't influence anything
   P <- matrix(kronecker(rowSums(adj),rowSums(adj))/(sum(cs)*sum(cs)),NROW(adj),NROW(adj))
   diag(P) <- 0
   
   # put the pieces together to make them more friendly
   B <- A-P

   # calculate the modularity based on being in the same group and the above contributions to modularity
   # we divide by two since each matrix is double the fun
   Q = (1/2.) * sum(Smat %*% t(Smat) * B)
   
   return(Q)
}

bLPW <- function (x,as.adjacency=FALSE) {
   # save the row order for later
   OrderVec <- rownames(x)

   # shuffle the adjacency so that its order doesn't influence anything
   x <- x[sample(c(1:NROW(x))),]
   x <- x[,sample(c(1:NCOL(x)))]

   # starting labels for the row communities
   lA <- c(1:NROW(x))
   names(lA) <- rownames(x)

   # build the Smat for modularity calculation
   Smat <- modulesToSmat(lA)

   ## Seeding the 'initial' modularity
   oldQ <- bQW(x,Smat)
   maxQ <- oldQ
   maxl <- lA
   
   # how's progress coming along?
   cat('Label propagation starting\n')
   cat('\rCurrent Q\t',oldQ)
   
   # Label propagation
   # # first, make sure that we get at least one increase in modularity
   Msteps <- 0
   while(TRUE){
      Msteps <- Msteps + 1

      lB <- lA

      ## update lB from lA in random order
      for(lsp in sample(rownames(x)))
      {
         # how strongly do other nodes' interactions overlap with those of lsp
         overlap <- x[lsp,] %*% t(x)

         # remove the node itself (maybe this isn't necessary?)
         # overlap <- overlap[,colnames(overlap) != lsp]

         # select the node with the greatest overlap
         donor <- names(overlap)[which.is.max(overlap)]

         # set the label of lsp to be that of the donor
         lB[lsp] <- lB[donor]
      }

      # build the new Smat for modularity calculation
      Smat <- modulesToSmat(lB)

      # what is the new modularity?
      newQ <- bQW(x,Smat)

      if(newQ > maxQ){
         maxQ <- newQ
         maxl <- lB
      }

      cat('\rCurrent Q\t',newQ)

      if(newQ > oldQ){
         break
      }
   }
   lA <- lB   

   # keep updating until there are no further changes
   Nsteps <- 1
   while(TRUE)
   {
      Nsteps <- Nsteps + 1
      lB <- lA

      ## update lB from lA in random order
      for(lsp in sample(rownames(x))){
         # how strongly do other nodes' interactions overlap with those of lsp
         overlap <- x[lsp,] %*% t(x)

         # remove the node itself
         # overlap <- overlap[,colnames(overlap) != lsp]

         # select the node with the greatest overlap
         donor <- names(overlap)[which.is.max(overlap)]

         # set the label of lsp to be that of the donor
         lB[lsp] <- lB[donor]
      }

      # build the new Smat for modularity calculation
      Smat <- modulesToSmat(lB)

      # what is the new modularity?
      newQ <- bQW(x,Smat)

      if(newQ > maxQ){
         maxQ <- newQ
         maxl <- lB
      }

      # progress is continuing
      cat('\rCurrent Q\t',newQ)

      # oops, we might not want to continue propagating!
      if(sum(lA != lB) == 0){
         break
      }else{
         lA <- lB
      }
   }
   
   # progress is coming along dandy!
   cat('\nLabel propagation ended after',Msteps,'+',Nsteps,'step(s)\n')
   
   # pretty it all up
   modules <- maxl
   modules <- as.numeric(as.factor(modules))
   names(modules) <- names(maxl)

   return(modules[OrderVec])
}

bBRIMW = function(x)
{
   # start with label propagation
   modules <- bLPW(x)

   # build the Smat for modularity calculation
   Smat <- modulesToSmat(modules)

   ## Some important values for modularity calculation
   # calcluate the column sums in the interaction matrix for use later
   cs <- colSums(x)
   
   # calculate the positive contribution to modularity of matching any two rows
   # also ignore the diagonal since it doesn't influence anything
   A <- (x %*% t(x)) / sum(cs*(cs-1))
   diag(A) <- 0
   
   # calculate the expected contribution to modularity that would occur if any two rows show overlapped by chance
   # also ignore the diagonal since it doesn't influence anything
   P <- matrix(kronecker(rowSums(x),rowSums(x))/(sum(cs)*sum(cs)),NROW(x),NROW(x))
   diag(P) <- 0
   
   # put the pieces together to make them more friendly
   B <- A-P

   # calculate the pre-BRIM modularity
   preBM <- (1/2.) * sum(Smat %*% t(Smat) * B)

   ## Optimization loop
   Nsteps <- 0
   cat('\nBRIM optimization starting\n')
   while(TRUE)
   {
      Nsteps <- Nsteps + 1

      # print it out for posterity's sake
      cat('\rCurrent Q\t',preBM)

      ## Product matrix for T & R
      rBT <- B %*% Smat

      ## Optimization
      cSmat <- Smat
      cSmat[,] <- 0
      cSmat[cbind(1:NROW(x),apply(rBT,1,which.max))] <- 1

      ## New bipartition
      cBM <- (1/2.) * sum(cSmat %*% t(cSmat) * B)

      if(cBM < preBM){
         break
      }else{
         Smat <- cSmat[,colSums(cSmat)>0]
         preBM <- cBM
      }
   }
   cat('\nBRIM convergence reached in',Nsteps,'step(s)\n')
   return(list(M=x,Smat=Smat,Q=preBM,c=NCOL(Smat)))
}

findModulesW = function(M,iter=50,cpu=1)
{
   usePar <- ifelse(cpu>1,TRUE,FALSE)
   sfInit(parallel=usePar,cpus=cpu,type='SOCK')
   if(usePar)
   {
      sfLibrary(nnet)
      sfExportAll()
   }
   if(is.null(rownames(M))) rownames(M) <- paste('r',c(1:NROW(M)),sep='')
   if(is.null(colnames(M))) colnames(M) <- paste('c',c(1:NCOL(M)),sep='')
   ModuleOutput <- sfLapply(c(1:iter),function(x) bBRIMW(M))
   sfStop()	
   return(ModuleOutput)
   #Qs <- unlist(lapply(ModulOutput,function(x)x$Q))
   #maxQs <- which.is.max(Qs)
   #return(ModulOutput[[maxQs]])
}
