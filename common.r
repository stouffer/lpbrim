
spread = function (v, m = 0, M = 1) 
{
   v <- v - min(v)
   v <- v/max(v)
   v <- v * (M - m)
   v <- v + m
   return(v)
}

mostFrequent = function(vec, w=NA)
{
   #if(is.na(w)) w <- rep(1, length(vec)) # what is this supposed to do? w isn't used anywhere else in the function!
   tvec <- table(vec)
   nvec <- as.vector(tvec[as.factor(vec)])
   return(which.is.max(nvec))
}

