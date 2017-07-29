#' A Plotting Function
#'
#' Produces a plot of the area-under-the-curve for the mortality data
#'
#' @param n the length of the vector Defaults to TRUE.
#' @param x the vector arguement.
#' @param add whether to add lines. Default is FALSE
#'
#' @keywords lplot
#'
#' @export
#'
#' @examples m1 <- Mortality$D.Male[which(Mortality$Year == 2008)]
#' @examples m2 <- Mortality$E.Male[which(Mortality$Year == 2008)]
#' @examples male.1 <- m1/m2
#' @examples male.2 <- log(male.1[!is.na(male.1)])
#' @examples lplot(1:length(male.2),male.2)
#'
#' @examples lplot(1:length(male.2),male.2,add=TRUE)
#'
#' @importFrom graphics plot abline hist lines par points text
#'
#' @importFrom stats density dnorm dt optim optimize pchisq pnorm
#' pt qexp qnorm qt rchisq rnorm rpois rt runif var



lplot=function(n,x,add=F){

	nn=length(n)
	ta=array(n,c(nn,3))
	xa=ta*0; xa[,2]=x
	ta=t(ta)
	xa=t(xa)
	dim(ta)=NULL
	dim(xa)=NULL
	if(add==F){ plot(ta,xa,type="l",lwd=0.5) }
	if(add==T){ lines(x,lwd=0.5,col=2) }
}

#' A Plotting Function
#'
#' Produces a plot of the area-under-the-curve for the mortality data, but lplot() inverted
#'
#'
#' @param n the length of the vector Defaults to TRUE.
#' @param x the vector arguement.
#' @param add whether to add lines. Default is FALSE
#'
#'
#' @keywords iplot
#'
#' @export
#'
#' @examples m1 <- Mortality$D.Male[which(Mortality$Year == 2008)]
#' @examples m2 <- Mortality$E.Male[which(Mortality$Year == 2008)]
#' @examples male.1 <- m1/m2
#' @examples male.2 <- log(male.1[!is.na(male.1)])
#' @examples iplot(1:length(male.2),male.2)
#' 
#' @examples iplot(1:length(male.2),male.2,add=TRUE)
#'
#' @importFrom graphics plot

iplot=function(n,x,add=F){

	nn=length(n)
	ta=array(n,c(nn,3))
	xa=ta*0; xa[,2]=x
	ta=t(ta)
	xa=t(xa)
	dim(ta)=NULL
	dim(xa)=NULL
	if(add==F){ plot(ta,-xa,type="l",lwd=0.5) }
	if(add==T){ lines(-x,lwd=0.5) }
}

#' A Plotting Function
#'
#' Produces a similar plot as lplot(), only a transposition of ages is made
#'
#'
#' @param n the length of the vector Defaults to TRUE.
#' @param x the vector arguement.
#' @param add whether to add lines. Default is FALSE
#'
#' @keywords vplot
#'
#' @export
#'
#'
#' @examples m1 <- Mortality$D.Male[which(Mortality$Year == 2008)]
#' @examples m2 <- Mortality$E.Male[which(Mortality$Year == 2008)]
#' @examples male.1 <- m1/m2
#' @examples male.2 <- log(male.1[!is.na(male.1)])
#' @examples vplot(1:length(male.2),male.2)
#'
#' @importFrom graphics plot


vplot=function(n,x,add=F){

	nn=length(n)
	ta=array(n,c(nn,3))
	xa=ta*0; xa[,2]=x
	ta=t(ta)
	xa=t(xa)
	dim(ta)=NULL
	dim(xa)=NULL
	if(add==F){ plot(-ta,-xa,type="l",lwd=0.5) }
	if(add==T){ lines(-ta,-xa,lwd=0.5,col=2) }
}

myplotsize=function(x0, x1, y0, y1, l = ""){
	plot(c(x0, x1), c(y0, y1), log = l, xlab = "", ylab = "", type = "n")
}


#' A Plotting Function
#'
#' Similar to pccopula(), but suitable when the dependence is
#' stronger at the older ages
#'
#'
#' @param theta gives the order.
#' @param pl gives the association, with a correction for the direction of dependence
#' @param z the length of the z axis Defaults to 10.
#' 
#' @keywords pgcopula()
#' @export
#' @examples #Examples
#'
#' @examples pgcopula(theta=1.3,pl=2,z=10)
#'
#'
#' @importFrom graphics plot persp contour title
#'


pgcopula=function(theta,pl=1,z){
	u1=((1:100)-0.5)/100
	u2=u1
	za=array(0,c(100,100))
	for(i  in 1:100){
		for(j in 1:100){
			v1=u1[i]
			v2=u2[j]
			x=(-log(v1))^theta+(-log(v2))^theta
			pd1=-theta/v1*((-log(v1))^(theta-1))
			pd2=-theta/v2*((-log(v2))^(theta-1))
			za[i,j]=pd1*pd2*(1/theta^2)*((theta-1)*(x^(1/theta-2))+x^(2/theta-2))*exp(-x^(1/theta))
		}
	}
	if(pl==1){
		myplotsize(0,1,0,1)
		vv=c((1:10)/10,1.5,2:10)
		contour(u1,u2,za,levels=vv,add=T)
		title(xlab="u_1",ylab="u_2")
	}
	else{
		nv=(0:33)*3+1
		persp(u1[nv],u1[nv],za[nv,nv],lwd=0.25,d=5,xlim=c(0,1),ylim=c(0,1),zlim=c(0,z),expand=0.75,xlab="u1",ylab="u2",zlab="copula density c(u1,u2)",ticktype="detailed",theta=300,phi=10)
	}
	a=paste("Gumbel copula density\nTheta = ",round(theta,2),sep="",collapse="")
	title(a)
}

#' A Plotting Function
#'
#' Produces a plot of a copula, which can be used to assess the dependency between
#' two sexes bounded by the actual and the expanded mortality estimates
#'
#' @param theta gives the order.
#' @param pl gives the association.
#' @param z the length of the z axis Defaults to 10.
#' 
#' @keywords pccopula()
#' @export
#' @examples #Examples
#'
#' @examples pccopula(theta=3,pl=.5,z=10)
#'
#'
#' @importFrom graphics plot persp contour title

	
pccopula=function(theta,pl=1,z){
	u1=((1:100)-0.5)/100
	u2=u1
	za=array(0,c(100,100))
	for(i  in 1:100){
		for(j in 1:100){
			za[i,j]=(theta+1)*((u1[i]*u2[j])^(-theta-1)) * ((u1[i]^(-theta)+u2[j]^(-theta)-1)^(-1/theta-2))
		}
	}
	if(pl==1){
		myplotsize(0,1,0,1)
		vv=c((1:10)/10,1.5,2:10)
		contour(u1,u2,za,levels=vv,add=T)
		title(xlab="u_1",ylab="u_2")
	}
	else{
		nv=(0:33)*3+1
		persp(u1[nv],u1[nv],za[nv,nv],lwd=0.25,d=5,xlim=c(0,1),ylim=c(0,1),zlim=c(0,z),expand=0.75,xlab="u1",ylab="u2",zlab="copula density c(u1,u2)",ticktype="detailed",theta=300,phi=10)
	}
	a=paste("Clayton copula density\nTheta = ",round(theta,2),sep="",collapse="")
	title(a)
}


#' A Plotting Function
#'
#' Produces a plot of the difference between the area-under-the-curve for the mortality data and the extended mortality boundaries
#'
#' @param n the length of the vector Defaults to TRUE.
#' @param x the vector arguement.
#' @param young the age at which the accident hump begins. Must be entered
#' @param old age at which, either mortality experience between males and females converge, or rapid acceleration of mortality. This is typically over 80 years.
#'
#' @keywords mmplot
#' @export
#' @examples #Examples
#' @examples m1 <- Mortality$D.Male[which(Mortality$Year == 2008)]
#' @examples m2 <- Mortality$E.Male[which(Mortality$Year == 2008)]
#' @examples male.1 <- m1/m2
#' @examples male.2 <- log(male.1[!is.na(male.1)])
#' @examples lplot(1:length(male.2),male.2)
#'
#'
#' @examples mmplot(1:length(male.2),male.2,young=17,old=80)
#'
#' @importFrom graphics plot


mmplot=function(n,x,young,old){

	nn=length(n)
	ta=array(n,c(nn,3))
	xa=ta*0; xa[,2]=x
	ta=t(ta)
	xa=t(xa)
	b=c(tan(tan(1)+((0.95-tan(1))/(young-old)*(x[is.finite(x)]-old))),1)
	tan=log((exp(x[!is.na(x)]))^b)
	dim(ta)=NULL
	dim(xa)=NULL
	plot(ta,xa,type="l",lwd=0.5)
	lines(tan,lwd=0.5,col=2)
	lines(x,lwd=0.5,col=2)
}

