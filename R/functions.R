# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

read <- function(name) {
  x <- read.csv(name)
  return(x)
}

mybin <- function(iter,n, p){
  # make a matrix to hold the samples
  #initially filled with NA's
  sam.mat=matrix(NA,nr=n,nc=iter, byrow=TRUE)
  #Make a vector to hold the number of successes in each trial
  succ=c()
  for( i in 1:iter){
    #Fill each column with a new sample
    sam.mat[,i]=sample(c(1,0),n,replace=TRUE, prob=c(p,1-p))
    #Calculate a statistic from the sample (this case it is the sum)
    succ[i]=sum(sam.mat[,i])
  }
  #Make a table of successes
  succ.tab=table(factor(succ,levels=0:n))
  #Make a barplot of the proportions
  barplot(succ.tab/(iter), col=rainbow(n+1), main="Binomial simulation", xlab="Number of successes")
  succ.tab/iter
}

myncurve = function(mu, sigma,a){
  curve(dnorm(x,mean=mu,sd=sigma), xlim = c(mu-3*sigma, mu + 3*sigma))
  xcurve=seq(-1000,a,length=1000)
  ycurve=dnorm(xcurve,mu,sigma)
  polygon(c(-1000,xcurve,a),c(0,ycurve,0),col="Red")
  prob=pnorm(a,mean=mu,sd=sigma)-pnorm(-1000,mean=mu,sd=sigma)
  prob=round(prob,4)
  text(a/2,0.5*0.1, paste("Area = ", prob, sep=""))

}
data = function(name){
  x <- read.csv("FIREDAM.csv")
  x = name
}

mycltu=function(n,iter,a=0,b=10){
  ## r-random sample from the uniform
  y=runif(n*iter,a,b)
  ## Place these numbers into a matrix
  ## The columns will correspond to the iteration and the rows will equal the sample size n
  data=matrix(y,nr=n,nc=iter,byrow=TRUE)
  ## apply the function mean to the columns (2) of the matrix
  ## these are placed in a vector w
  w=apply(data,2,mean)
  ## We will make a histogram of the values in w
  ## How high should we make y axis?
  ## All the values used to make a histogram are placed in param (nothing is plotted yet)
  param=hist(w,plot=FALSE)
  ## Since the histogram will be a density plot we will find the max density

  ymax=max(param$density)
  ## To be on the safe side we will add 10% more to this
  ymax=1.1*ymax
  ## Now we can make the histogram
  hist(w,freq=FALSE,  ylim=c(0,ymax), main=paste("Histogram of sample mean",
                                                 "\n", "sample size= ",n,sep=""),xlab="Sample mean")
  ## add a density curve made from the sample distribution
  lines(density(w),col="Blue",lwd=3) # add a density plot
  ## Add a theoretical normal curve
  curve(dnorm(x,mean=(a+b)/2,sd=(b-a)/(sqrt(12*n))),add=TRUE,col="Red",lty=2,lwd=3) # add a theoretical curve
  ## Add the density from which the samples were taken
  curve(dunif(x,a,b),add=TRUE,lwd=4)

}

myboot2<-function(iter=10000,x,fun="mean",alpha=0.05,cx=1.5,...){  #Notice where the ... is repeated in the code
  n=length(x)   #sample size

  y=sample(x,n*iter,replace=TRUE)
  rs.mat=matrix(y,nr=n,nc=iter,byrow=TRUE)
  xstat=apply(rs.mat,2,fun) # xstat is a vector and will have iter values in it
  ci=quantile(xstat,c(alpha/2,1-alpha/2))# Nice way to form a confidence interval
  # A histogram follows
  # The object para will contain the parameters used to make the histogram
  para=hist(xstat,freq=FALSE,las=1,
            main=paste("Histogram of Bootstrap sample statistics","\n","alpha=",alpha," iter=",iter,sep=""),
            ...)

  #mat will be a matrix that contains the data, this is done so that I can use apply()
  mat=matrix(x,nr=length(x),nc=1,byrow=TRUE)

  #pte is the point estimate
  #This uses whatever fun is
  pte=apply(mat,2,fun)
  abline(v=pte,lwd=3,col="Black")# Vertical line
  segments(ci[1],0,ci[2],0,lwd=4)      #Make the segment for the ci
  text(ci[1],0,paste("(",round(ci[1],2),sep=""),col="Red",cex=cx)
  text(ci[2],0,paste(round(ci[2],2),")",sep=""),col="Red",cex=cx)

  # plot the point estimate 1/2 way up the density
  text(pte,max(para$density)/2,round(pte,2),cex=cx)

  invisible(list(ci=ci,fun=fun,x=x))# Some output to use if necessary
}

mymld=function(lfun,x,param,delta=0.0001){    # param = parameter values, delta=accuracy, x=data
  z=outer(x,param,lfun)    # create outer product and evaluate at lfun
  y=apply(z,2,sum) # x by param, 2=columns , sum columns = sum of log lik

  i=max(which(y==max(y)))# the index for which y is biggest, if two then take the last one

  param2=seq(param[i-2],param[i+2],by=delta)# The maximum will be between these two, increments by delta
  zz=outer(x,param2,lfun) # new z, call it zz
  yy=apply(zz,2,sum)   # new y, call it yy
  ii=max(which(yy==max(yy)))# new i,  call it ii , if two, take max of them (last one)
  layout(matrix(c(1,2),nr=1,nc=2,byrow=TRUE))# divide plotting space for two graphs
  plot(param,y,col="Blue",type="l",lwd=2,ylab="Log. Lik.",xlab=expression(theta))# plot log lik Vs parameter values
  abline(v=param[i],lwd=2,col="Red") # Show vertical line at estimated value
  axis(3,param[i],round(param[i],2))
  points(param[i],y[i],pch=19,cex=1.5,col="Black")# Plot the point
  plot(param2,yy,col="Blue",type="l",lwd=2,ylab="Log. Lik.",xlab=expression(theta),las=2) # construct new plot for refined estimate
  abline(v=param2[ii],lwd=2,col="Red")  # new verical line
  val=round(param2[ii],abs(log10(delta))) ## rounds to the nth place where n is st delta=10^-n.
  axis(3,param2[ii],val)
  points(param2[ii],yy[ii],pch=19,cex=1.5,col="Black")
}
