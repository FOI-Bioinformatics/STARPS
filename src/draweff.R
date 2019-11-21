draweff<-function(size=1.8E6,p1=0.4,p2=0.5,mu=0.1){
  # draw effect of a mutation on fitness
  # p1 p2
  # mu is the mean of the exponential distribution or here:
  # the average selective advantage (or disadvantage) of mutations
  # in Rozen et al. 2002, mu = 1/35 was used
  eff<-vector(size,mode = "numeric")
  cdf<-c(p1,p1+p2,1) # define limits for generating sign of fitness effect
  
  # step 1: decide which distribution to use 

  ra<-1/mu
  u <- runif(size,0.0,1.0)
  # step 2: make a draw from that distribution and repeat it size times
  for(i in 1:size){
 #   u <- runif(1,0.0,1.0)
   
    if(u[i] <= cdf[1])
      eff[i] <- - rexp(1,rate=ra) 
    else if(u[i] > cdf[2])
      eff[i] <- rexp(1,rate=ra) 
  }
  eff

}

