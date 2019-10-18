#Returns the dpc distance for persistence disgrams d1 and d2 for Betti number beta

infnorm <- function(x,y){
  return(max(abs(x[1]-y[1]),abs(x[2]-y[2])))
}

schudist <- function(d1,d2, beta,c,p){
  library('clue')
  d1 = as.matrix(d1)
  d2 = as.matrix(d2)
  
  wd1 = d1[which(d1[,1]==beta),2:3, drop = FALSE]
  wd2 = d2[which(d2[,1]==beta),2:3, drop = FALSE]
  n = dim(wd1)[1]
  m = dim(wd2)[1]
  
  
  if(n == 0 || m == 0)
  {
    if(n == 0 && m == 0){
      return(0)
    }
    return(c)
  }
  
  
  #always make wd1 smaller
  if(m < n){
    dum1 = n
    dum2 = wd1
    n = m
    m = dum1
    wd1 = wd2
    wd2 = dum2
  }
  #Create distance matrix
  distmat = matrix(nrow=n,ncol=m)
  for(i in 1:n){
    for(j in 1:m){
      distmat[i,j] = min(c,infnorm(wd1[i,],wd2[j,]))^p
    }
  }
  solution = solve_LSAP(distmat)
  cost = 0
  for(k in 1:length(solution)){
    cost = cost + distmat[k,solution[k]]
  }
  cost = cost + (c^p)*(m-n)
  cost = cost*(1/m)
  cost = cost^(1/p)
  return(cost)
}

