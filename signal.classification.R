#install.packages("TDA")
#install.packages("ggplot2")
#install.packages("clue")
#install.packages("caret")
library(clue)
library(TDA)
library(ggplot2)
library(caret)

source("~/Downloads/dpc.R")

set.seed(40)
tau <- 45
p <- 3
c <- 2

#Classification of bistable and random walk signals using dpc
#distance

LinSigGaussNoise <- function(initialx, noisemode){
  # Generates a linear signal with Gaussian noise
  # Noise distribution is N(noisemode, noisestd)
  # Initial point for signal is (0, initialx,initialy)
  # Input: initial-float , noisemode-float, noisevar-float
  # Output: signal- 2x200 matrix
  time     <- 0:200
  amps     <- integer(201)
  amps[1] <- initialx
  noisestd <- runif(1,0,0.05)
  for (i in 1:200){
    amps[i+1] = amps[i]+ rnorm(1,mean = noisemode,sd = noisestd)
    end
  }
  
  rwsignal = cbind(time,amps)
  return(rwsignal)
}

Bistable <- function(initialx,alpha){
  # Generates bistable signal acc to Andy's paper p14
  # Initial point for signal is (0, initialx,initialy)
  # Input: initial , alpha, noisestd-float
  # Output: signal- 2x200 matrix
  time     <- 0:200
  amps     <- integer(201)
  amps[1] <- initialx
  noisestd <- runif(1,0,0.05)
  for (i in 1:200){
    amps[i+1] = alpha*sin(amps[i]) + rnorm(1,mean = 0,sd = noisestd)
    end
  }
  
  bisignal = cbind(time,amps)
  return(bisignal)
}



#Generate two sample signals
rwsignal <- LinSigGaussNoise(0,0)
bisignal <- Bistable(0,2.5)

par(mfrow = c(1,2))
plot(rwsignal)
plot(bisignal)

#Look at autocorrelation for first zero
acf(rwsignal[,2],lag.max = 100)
acf(bisignal[,2],lag.max = 100)
tau <- 2

twodPointCloud <- function(signal,tau){
  # Computes the 3-dimensionial point cloud for a signal
  # Input: 2x201 matrix
  # Output: 3x200 matrix
  
  amps <- signal[,2]
  end <- 201-tau
  begin <- 1 + tau
  x1 <- amps[1:end]
  x2 <- amps[begin:201]
  pointcloud <- cbind(x1,x2)

}

#Make sin signal and plot point cloud
time     <- 0:200
amps     <- integer(201)
amps <- sin(time)
sinsignal = cbind(time,amps,amps)
pointcloud <- twodPointCloud(sinsignal,1)
par(mfrow = c(1,1))
plot(pointcloud) #ellipse


#Create point clouds and persistence diagrams
rwpointcloud <- twodPointCloud(rwsignal,tau)
rw_diagram <- ripsDiag(rwpointcloud,1,5)
bipointcloud <- twodPointCloud(bisignal,tau)
bi_diagram <- ripsDiag(bipointcloud,1,5)
par(mfrow = c(1,2))
plot(rwpointcloud)
plot(bipointcloud)


#plot signals and the persistance diagrams
par(mfrow = c(2, 2))
plot(rwsignal)
plot(bisignal)
plot(rw_diagram[["diagram"]])
plot(bi_diagram[["diagram"]])



Dis_to_class <- function(diag,class,dim,trainsize){
  sum <- 0
  for(i in 1:trainsize){
    if(class == 1){
      first_diag <- rwTraining[[i]]
      sum = sum + schudist(first_diag[["diagram"]],diag[["diagram"]]
                           ,dim,p,c)
    }
    if(class == 2){
      first_diag <- biTraining[[i]]
      sum = sum + schudist(first_diag[["diagram"]],diag[["diagram"]]
                           ,dim,p,c)
    }
  }
  average_dist <- sum/trainsize
  return(average_dist)
}

Classify <- function(diag,trainsize){
  Dis <- matrix(0,2,2)
  bsum <- integer(2)
  for(i in 1:2){
    Dis[1,i] <- Dis_to_class(diag,i,0,trainsize)
    Dis[2,i] <- Dis_to_class(diag,i,1,trainsize)
  }
  bsum[1] <- sum(Dis[,1])
  bsum[2] <- sum(Dis[,2])
  label <- which(bsum == min(bsum))
  return(label)
}

#Make sets 
ss <- 10 #make div by 10
rwData <- integer(ss)
biData <- integer(ss)
bothdata = list()
for(i in 1:ss){
  rwsignal <- LinSigGaussNoise(0,0) 
  bisignal <- Bistable(0,2.5) 
  rwpointcloud <- twodPointCloud(rwsignal,tau)
  bipointcloud <- twodPointCloud(bisignal,tau)
  bothdata[[i]] <- ripsDiag(rwpointcloud,1,5)
  bothdata[[i+ss]] <- ripsDiag(bipointcloud,1,5)
}

a <- rep(1,ss)
b <- rep(2,ss) 
actuallabels <- as.factor(c(a,b))

folds = createFolds(actuallabels,k = 10)

testlabels = integer(2*ss)
for(i in 1:10){
  Testing  <- bothdata[folds[[i]]]
  rwTraining = bothdata[c(1:ss)[-folds[[i]]]]
  biTraining = bothdata[c((ss+1):(2*ss))[-folds[[i]]]]
  testlabels[folds[[i]]] <- lapply(Testing,trainsize = ss-(ss/10),Classify)
}

testlabels = as.factor(unlist(testlabels))
confusion <- confusionMatrix(data = testlabels,actuallabels)
confusion



