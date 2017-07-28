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

par1=function(){
	par(mfrow = c(1, 1))
	par(mar = c(6, 4, 2, 1), lwd = 1, xaxs = "i", yaxs = "i", tck = -0.02, err = -1.)
	0
}


setaxes=function(x0, x1, y0, y1, l = ""){
	plot(c(x0, x1), c(y0, y1), log = l, xlab = "", ylab = "", type = "n")
}

cdfplot=function(x,a=0,co=1,lt=1){

	x=x[order(x)]	
	n=length(x)
	nv=(1:n)*2
	xv=nv*0
	xv[nv]=x
	xv[nv-1]=x
	yv=xv*0
	yv[nv]=(1:n)/n
	yv[nv-1]=(0:(n-1))/n
	if(a==0){ setaxes(x[1],x[n],0,1) }
	lines(xv,yv,col=co,lty=lt)
}

	
pdfplot=function(x,a=0,co=1,lt=1){

	res=density(x)
	if(a==0){ 
		setaxes(min(res$x),max(res$x),0,max(res$y)*1.2)
	}
	lines(res$x,res$y,col=co,lty=lt)
}


chi2test.normal=function(d,nbin=100,pr=0){

	n=length(d)

	mu=mean(d)
	sd=sqrt(mean((d-mu)^2))
	uv=(0:nbin)/nbin

	br=qnorm(uv,mu,sd)
	obs=(1:nbin)*0
	for(e in 1:nbin){

		obs[e]=sum((d > br[e])*(d <= br[e+1]))
	}
	ex=n/nbin	
	cs=sum((obs-ex)^2/ex)
	df=nbin-3
	pv=1-pchisq(cs,df)
	if(pr == 1){
		uh=pnorm(d,mu,sd)
		hist(uh,breaks=uv)
		abline(h=ex,col=2)
	}
	cat("Chi-squared statistic is ",cs, "with ",df," degrees of freedom.\n")
	cat("The p-value for this is ",pv,".\n")
}



jbtest=function(x){

	m1=mean(x)
	n=length(x)
	m2=var(x)*(n-1)/n	
	m3=mean((x-m1)^3)/(m2^1.5)
	m4=mean((x-m1)^4)/(m2^2)

	jb=1/6*n*(m3*m3+0.25*(m4-3)^2)

	pv=1-pchisq(jb,df=2)
	cat("   Skewness =",m3,"\n")
	cat("   Kurtosis =",m4,"\n")
	cat("Jarque-Bera =",jb,"\n")
	cat("    p-value =",pv,"\n")
	pv
}




kurtosis=function(x){

	m1=mean(x)
	n=length(x)
	m2=var(x)*(n-1)/n	
	m3=mean((x-m1)^3)/(m2^1.5)
	m4=mean((x-m1)^4)/(m2^2)
	cat("   Kurtosis =",m4,"\n")
	m4
}

	

skewness=function(x){

	m1=mean(x)
	n=length(x)
	m2=var(x)*(n-1)/n	
	m3=mean((x-m1)^3)/(m2^1.5)
	m4=mean((x-m1)^4)/(m2^2)
	cat("   Skewness =",m3,"\n")
	m3
}





rnct=function(n,df,ncp=0){

	if(df > 2){ 
		mu=sqrt(df/2)*gamma(0.5*(df-1))/gamma(0.5*df)*ncp
		mu2=df/(df-2)*(1+ncp*ncp)
		sd=sqrt(mu2-mu^2)
		z=rt(n,df=df,ncp=ncp)
		Z=(z-mu)/sd
		Z
	} else { cat("Error: parameter df must be > 2\n") }
}


dnct=function(x,df,ncp=0){

	mu=sqrt(df/2)*gamma(0.5*(df-1))/gamma(0.5*df)*ncp
	mu2=df/(df-2)*(1+ncp*ncp)
	sd=sqrt(mu2-mu^2)
	density=sd*dt(mu+sd*x,df=df,ncp=ncp)
	density
}

pnct=function(x,df,ncp=0){

	mu=sqrt(df/2)*gamma(0.5*(df-1))/gamma(0.5*df)*ncp
	mu2=df/(df-2)*(1+ncp*ncp)
	sd=sqrt(mu2-mu^2)
	prob=pt(mu+sd*x,df=df,ncp=ncp)
	prob
}

qnct=function(u,df,ncp=0){

	mu=sqrt(df/2)*gamma(0.5*(df-1))/gamma(0.5*df)*ncp
	mu2=df/(df-2)*(1+ncp*ncp)
	sd=sqrt(mu2-mu^2)
	quant=(qt(u,df=df,ncp=ncp)-mu)/sd
	quant
}



ll.garch.normal=function(p,x,summary=0,print=1){

	sigma1=exp(p[1])
	alpha0=exp(p[2])
	alpha1=exp(p[3])
	beta=exp(p[4])
	mu=p[5]
	n=length(x)
	vol=(1:n)*0
	vol[1]=sigma1^2
	for(i in 2:n){ vol[i]=alpha0+alpha1*(x[i-1]-mu)^2+beta*vol[i-1] }
	ll=sum(-0.5*log(vol)+dnorm((x-mu)/sqrt(vol),log=T))
	sigma=sqrt(vol)
	z=(x-mu)/sigma
	if(print==1){ cat(ll,"\n") }
	if(summary == 0){ ll } else { list(sigma=sigma,Z=z,ll=ll,x=x) }
}


ll.garch.t=function(p,x,summary=0,print=1){

	sigma1=exp(p[1])
	alpha0=exp(p[2])
	alpha1=exp(p[3])
	beta=exp(p[4])
	mu=p[5]
	df=p[6]
	n=length(x)
	vol=(1:n)*0
	vol[1]=sigma1^2
	for(i in 2:n){ vol[i]=alpha0+alpha1*(x[i-1]-mu)^2+beta*vol[i-1] }
	ll=sum(-0.5*log(vol)+log(dnct((x-mu)/sqrt(vol),df=df)))
	sigma=sqrt(vol)
	z=(x-mu)/sigma
	if(print==1){ cat(ll,"\n") }
	if(summary == 0){ ll } else { list(sigma=sigma,Z=z,ll=ll,x=x) }
}



ll.garch.nct=function(p,x,summary=0,print=1){

	sigma1=exp(p[1])
	alpha0=exp(p[2])
	alpha1=exp(p[3])
	beta=exp(p[4])
	mu=p[5]
	df=p[6]
	ncp=p[7]
	n=length(x)
	vol=(1:n)*0
	vol[1]=sigma1^2
	for(i in 2:n){ vol[i]=alpha0+alpha1*(x[i-1]-mu)^2+beta*vol[i-1] }
	ll=sum(-0.5*log(vol)+log(dnct((x-mu)/sqrt(vol),df=df,ncp=ncp)))
	sigma=sqrt(vol)
	z=(x-mu)/sigma
	if(print==1){ cat(ll,"\n") }
	if(summary == 0){ ll } else { list(sigma=sigma,Z=z,ll=ll,x=x) }
}


fit.garch11.normal=function(x,p0=-1,print=1){
	if(p0==-1){ p0=c(0.01,0.000001,0.09,0.9,0) }
	p0[1:4]=log(p0[1:4])
	res=optim(p0,ll.garch.normal,x=x,print=print,control=list(fnscale=-1,maxit=10000))	# (***)

	res1=ll.garch.normal(res$par,x,summary=1,print=print)
	param=res$par	
	param[1:4]=exp(param[1:4])	
	names(param)=c("sigma(1)","alpha0","alpha1","beta1","mu")
	res1$par=param
	res1
}


fit.garch11.t=function(x,p0=-1,print=1){
	if(p0[1] == -1){ p0=c(0.01,0.000001,0.09,0.9,0,8) }
	p0[1:4]=log(p0[1:4])
	res=optim(p0,ll.garch.t,x=x,print=print,control=list(fnscale=-1,maxit=10000))	# (***)
	res1=ll.garch.t(res$par,x,summary=1,print=print)
	param=res$par	
	param[1:4]=exp(param[1:4])	
	names(param)=c("sigma(1)","alpha0","alpha1","beta1","mu","DF")
	res1$par=param
	res1
}


fit.garch11.nct=function(x,p0=-1,print=1){
	if(p0 == -1){ p0=c(0.01,0.000001,0.09,0.9,0,8,0) }
	if(p0 == -2){ res=fit.garch11.t(x,print=print); p0=res$par; attributes(p0)=NULL; p0=c(p0,0) }
	p0[1:4]=log(p0[1:4])
	res=optim(p0,ll.garch.nct,x=x,print=print,control=list(fnscale=-1,maxit=10000))	# (***)
	res1=ll.garch.nct(res$par,x,summary=1,print=print)
	param=res$par	
	param[1:4]=exp(param[1:4])	
	names(param)=c("sigma(1)","alpha0","alpha1","beta1","mu","DF","ncp")
	res1$par=param
	res1
}


update.sigma=function(sigma0,alpha0,alpha1,beta1,Z0){
	sigma1=sqrt(alpha0+(alpha1*Z0^2+beta1)*sigma0^2)
	sigma1
}


update.sigma.R=function(res){
	Z=res$Z
	n=length(Z)
	Z0=Z[n]
	sigma0=res$sigma[n]
	param=res$par
	attributes(param)=NULL
	alpha0=param[2]
	alpha1=param[3]
	beta1=param[4]
	sigma1=sqrt(alpha0+(alpha1*Z0^2+beta1)*sigma0^2)
	sigma1
}


simulate.garch11.normal=function(n.sim,X0,T1,sigma1,alpha0,alpha1,beta1,mu){
	X.array=array(0,c(n.sim,T1+1))	
	d.array=array(0,c(n.sim,T1))	
	V.array=array(0,c(n.sim,T1))	
	Z.array=array(rnorm(n.sim*T1),c(n.sim,T1))
	V.array[,1]=sigma1^2
	for(i in 2:T1){
		V.array[,i]=alpha0+(alpha1*Z.array[,i-1]^2+beta1)*V.array[,i-1]
	}
	sigma.array=sqrt(V.array)
	X.array[,1]=X0
	for(i in 1:T1){
		d.array[,i]=mu+sigma.array[,i]*Z.array[,i]
		X.array[,i+1]=X.array[,i]*exp(d.array[,i])
	}
	list(t=0:T1,X=X.array,d=d.array,sigma=sigma.array,Z=Z.array)
}


simulate.garch11.normal.R=function(n.sim,X0,T1,res){
	sigma1=update.sigma.R(res)
	param=res$par
	attributes(param)=NULL
	alpha0=param[2]
	alpha1=param[3]
	beta1=param[4]
	mu=param[5]
	Xres=simulate.garch11.normal(n.sim,X0,T1,sigma1,alpha0,alpha1,beta1,mu)
	Xres
}


qqt=function(xv,pl=1,mle=0){

	n=length(xv)
	m1=mean(xv)
	xv=xv[order(xv)]
	xv1=xv-m1
	m2=var(xv1)*(n-1)/n
	m4=mean(xv1^4)
	m5=m4/m2/m2/3
	df=(4*m5-2)/(m5-1)
	sd=sqrt(m2*(df-2)/df)
	if(mle==1){ pv=mle.t(xv); df=pv[3]; sd=pv[2]; m1=pv[1] }

	uv=((1:n)-0.5)/n
	xv3=qt(uv,df=df)
	if(pl == 1){
		plot(xv3,xv,xlab="",ylab="")
		title(xlab="Theoretical quantiles",ylab="Sample quantiles")
	}
	c(m1,sd,df)
}


mle.normal=function(x,pr=1){
        mu=mean(x)
        sigma=mean((x-mu)^2)^0.5
        ll=sum(log(dnorm(x,mu,sigma)))
        if(pr == 1){
                cat("mu    = ",mu,"\n")
                cat("sigma = ",sigma,"\n")
                cat("log-likelihood = ",ll,"\n")
        }
        c(mu,sigma)
}


mle.t=function(x,pr=1){
	options(warn=-1)	
	pv0=qqt(x,pl=0)
	if(pv0[3] <= 0){ pv0[3]=10 }
	pv=optim(pv0,ll.t,x=x,control=list(fnscale=-1))$par
	if(pr==1){
		cat("mu    =",pv[1],"\n")
		cat("sigma =",pv[2],"\n")
		cat("nu    =",pv[3],"\n")
		cat("log-likelihood =",ll.t(pv,x),"\n")
	}
	options(warn=0)
	pv
}


ll.t=function(p,x){
	n=length(x)
	mu=p[1]
	sigma=p[2]
	nu=p[3]
	y=(x-mu)/sigma
	llv=dt(y,df=nu,log=TRUE)-log(sigma)
	sum(llv)
}


qqnct=function(x){
        pv=mle.nct(x)
        mu=pv[1]
        sigma=pv[2]
        nu=pv[3]
        ncp=pv[4]
        x=x[order(x)]
        n=length(x)
        uv=((1:n)-0.5)/n
        y=x
        for(e in 1:10){
                y=y-(pt((y-mu)/sigma,df=nu,ncp=ncp)-uv)/(dt((y-mu)/sigma,df=nu,ncp=ncp)/sigma)
        }
        plot(y,x)
        pv
}


mle.nct=function(x,pr=1){
	pv0=c(mle.t(x,pr=0),0)
	options(warn=-1)	
	pv=optim(pv0,ll.nct,x=x,control=list(fnscale=-1))$par
	if(pr==1){
		cat("mu    =",pv[1],"\n")
		cat("sigma =",pv[2],"\n")
		cat("nu    =",pv[3],"\n")
		cat("ncp   =",pv[4],"\n")
		cat("log-likelihood =",ll.nct(pv,x),"\n")
	}
	options(warn=0)	
	pv
}

ll.nct=function(p,x){
	n=length(x)
	mu=p[1]
	sigma=p[2]
	nu=p[3]
	ncp=p[4]
	y=(x-mu)/sigma
	llv=dt(y,df=nu,ncp=ncp,log=TRUE)-log(sigma)
	sum(llv)
}


chi2test.t=function(d,nbin=100,pr=0){
	n=length(d)
	params=mle.t(d)
	mu=params[1]
	sigma=params[2]
	nu=params[3]
	uv=(0:nbin)/nbin
	br=qt(uv,nu)*sigma+mu
	obs=(1:nbin)*0
	for(e in 1:nbin){

		obs[e]=sum((d > br[e])*(d <= br[e+1]))
	}
	ex=n/nbin	
	cs=sum((obs-ex)^2/ex)	
	df=nbin-3-1	
	pv=1-pchisq(cs,df)
	if(pr == 1){
		uh=pt((d-mu)/sigma,nu)
		hist(uh,breaks=uv)
		abline(h=ex,col=2)
	}
	cat("Chi-squared statistic is ",cs, "with ",df," degrees of freedom.\n")
	cat("The p-value for this is ",pv,".\n")
}



chi2test.nct=function(d,nbin=100,pr=0){
	n=length(d)
	params=mle.nct(d)
	mu=params[1]
	sigma=params[2]
	nu=params[3]
	ncp=params[4]
	uv=pt((d-mu)/sigma,df=nu,ncp=ncp)

	obs=(1:nbin)*0
	for(e in 1:nbin){

		obs[e]=sum((uv > (e-1)/nbin )*(uv <= e/nbin ))
	}
	ex=n/nbin	
	cs=sum((obs-ex)^2/ex)
	df=nbin-4-1
	pv=1-pchisq(cs,df)
	if(pr == 1){
		hist(uv,breaks=(0:nbin)/nbin)
		abline(h=ex,col=2)
	}
	cat("Chi-squared statistic is ",cs, "with ",df," degrees of freedom.\n")
	cat("The p-value for this is ",pv,".\n")
}

qqpareto=function(xv,pl=1,mle=0){
	n=length(xv)
	xv=xv[order(xv)]
	m1=mean(xv)
	m2=mean(xv^2)
	alpha=2*(m1*m1-m2)/(2*m1*m1-m2)
	lambda=m1*(alpha-1)
	if(mle==1){ pv=mle.pareto(xv); alpha=pv[1]; lambda=pv[2] }
	uv=((1:n)-0.5)/n
	xv3=qpareto(uv,lambda,alpha)
	if(pl == 1){
		plot(xv3,xv,xlab="",ylab="")
		title(xlab="Theoretical quantiles",ylab="Sample quantiles")
	}
	c(alpha,lambda,-ll.pareto(c(alpha,lambda),xv))
}


mle.pareto=function(x,pr=1){
	pv0=qqpareto(x,pl=0,mle=0)[1:2]
#	pv0=c(4,1)
	pv=optim(pv0,ll.pareto,x=x)$par
	if(pr==1){
		cat("alpha  =",pv[1],"\n")
		cat("lambda =",pv[2],"\n")
	}
	pv
}


ll.pareto=function(p,x){
	alpha=p[1]
	lambda=p[2]
	n=length(x)
	llv=log(dpareto(x,lambda,alpha))
	-sum(llv)
}


ppareto=function(x,lambda,alpha){
	1-(lambda/(lambda+x))^alpha
}


qpareto=function(u,lambda,alpha){
	x=lambda/((1-u)^(1/alpha))-lambda
	x
}

rpareto=function(n,lambda,alpha){
	uv=runif(n)
	qpareto(uv,lambda,alpha)
}

dpareto=function(x,lambda,alpha){
	alpha*(lambda^alpha)/((lambda+x)^(alpha+1))
}

qqgev=function(x,xi){
	n=length(x)
	uv=((1:n)-0.5)/n
	if(xi > 0){
		qv=1/xi*((-log(uv))^(-xi)-1)
	}
	if(xi == 0){
		qv=-log(-log(uv))
	}
	if(xi < 0){
		qv=1/xi*((-log(uv))^(-xi)-1)
	}
	x=x[order(x)]
	x0=min(x)
	x1=max(x)
	q0=min(qv)
	q1=max(qv)
	setaxes(1.05*q0-0.05*q1,-0.05*q0+1.05*q1,1.05*x0-0.05*x1,-0.05*x0+1.05*x1)
	points(qv,x)
	title(xlab="Theoretical GEV quantiles",ylab="Observed quantiles")
}


ll.gev=function(pv,x){
	xi=pv[1]
	mu=pv[2]
	sigma=pv[3]
	if(xi > 0){
		y=(1+xi*(x-mu)/sigma)
		log.density=-(y^(-1/xi))-(xi+1)/xi*log(y)-log(sigma)
	}
	if(xi == 0){
		y=(x-mu)/sigma
		log.density=-exp(-y)-y-log(sigma)
	}
	if(xi < 0){
		y=(1+xi*(x-mu)/sigma)
		log.density=-(y^(-1/xi))-(xi+1)/xi*log(y)-log(sigma)
	}
	log.likelihood=sum(log.density)
	log.likelihood
}


mle.gev=function(x,pv0=c(0.1,0,0.1),pr=1){
	xi=pv0[1]
	sigma=pv0[3]
	if(xi > 0){ pv0[2]=min(x)+sigma/xi-1 }
	if(xi < 0){ pv0[2]=max(x)+sigma/xi+1 }
        pv=optim(pv0,ll.gev,x=x,control=list(fnscale=-1))$par
        if(pr==1){
                cat("xi  =",pv[1],"\n")
                cat("mu =",pv[2],"\n")
                cat("sigma =",pv[3],"\n")
        }
        pv
}



meplot=function(x,umin=0,a=0,co=1,ci=0){
	x=x[order(-x)]
	l=length(x)
	n0=sum(x > umin)
	umax=max(x)
	me=(1:n0)*0
	se=me
	confidence.interval=array(0,c(2,n0))
	l1=length(me)
	for(e in n0:2){
		y=x[1:(e-1)]-x[e]
		me[e]=sum(y)/(e-1)
		se[e]=sqrt(var(y)/(e-1))
		confidence.interval[,e]=me[e]+c(-1.96,1.96)*se[e]
	}
	me[1]=0	
	xv=x[n0:1]
	me=me[n0:1]
	x0=min(xv)
	x1=max(xv)
	if(a==0){ setaxes(x0,x1,0,max(me)*1.1) }
	points(xv,me,pch=16,col=co,cex=0.5)
	if(ci == 1){
		lines(xv,confidence.interval[1,n0:1],lwd=0.5)
		lines(xv,confidence.interval[2,n0:1],lwd=0.5)
	}
	x0=xv[1]
	y0=me[1]
	alv=c(2,3,4,5)
	ylocation=c(0.8,0.7,0.6,0.5)
	if(a==0){
		for(e in 1:4){
			abline(y0-x0*1/(alv[e]-1),1/(alv[e]-1),col=2,lwd=2)
			lab=paste("xi =1/",alv[e],sep="",collapse="")
			lscale=max(me)-y0
			a0=y0
			b0=1/(alv[e]-1)
			wy=y0+ylocation[e]*lscale
			wx=(wy-a0)/b0+x0
			text(wx,wy,lab)
		}
	}

	title(xlab="Threshold, u",ylab="Sample mean excess")
}



hill.plot=function(x,umin=0,lambda=0){
	x=x[order(-x)]
	n=sum(x > umin)
	y=log(lambda+x[1:n])
	alpha.hat=(1:n)*0
	xi.hat=(1:n)*0
	xi.se=(1:n)*0
	xi.ci=array(0,c(2,n))
	alpha.ci=array(0,c(2,n))
	n0=20
	for(e in n0:n){
		alpha.hat[e]=(mean(y[1:e])-y[e])^(-1)
		xi.hat[e]=(mean(y[1:e])-y[e])
		xi.se[e]=sqrt(1/e*var(y[1:e]))
		xi.ci[,e]=xi.hat[e]+c(-2,2)*xi.se[e]
		if(xi.ci[1,e] < 0) { xi.ci[1,e]=0.001 }
#		if(xi.ci[2,e] < 0) { xi.ci[1,e]=0.001 }
		alpha.ci[,e]=1/xi.ci[,e]
	}
	ra=max(alpha.ci)-min(alpha.ci)
	setaxes(0,n,min(alpha.ci)-0.05*ra,max(alpha.ci)+0.05*ra)
	lines(n0:n,alpha.hat[n0:n])
	lines(n0:n,alpha.ci[1,n0:n],lty=2)
	lines(n0:n,alpha.ci[2,n0:n],lty=2)
}

hill.est=function(x,umin=0,lambda=0){
	x=x[order(-x)]
	n=sum(x > umin)
	y=log(lambda+x[1:n])
	xi=mean(y)-y[n]
	alpha=1/xi
	list(xi=xi,alpha=alpha)
}


meplot=function(x,umin=0,a=0,co=1,ci=0){
	x=x[order(-x)]
	l=length(x)
	n0=sum(x > umin)
	umax=max(x)
	me=(1:n0)*0
	se=me
	confidence.interval=array(0,c(2,n0))
	l1=length(me)
	for(e in n0:2){
		y=x[1:(e-1)]-x[e]
		me[e]=sum(y)/(e-1)
		se[e]=sqrt(var(y)/(e-1))
		confidence.interval[,e]=me[e]+c(-1.96,1.96)*se[e]
	}
	me[1]=0		
	xv=x[n0:1]
	me=me[n0:1]
	x0=min(xv)
	x1=max(xv)
	if(a==0){ setaxes(x0,x1,0,max(me)*1.1) }
	points(xv,me,pch=16,col=co,cex=0.5)
	if(ci == 1){
		lines(xv,confidence.interval[1,n0:1],lwd=0.5)
		lines(xv,confidence.interval[2,n0:1],lwd=0.5)
	}
	x0=xv[1]
	y0=me[1]
	alv=c(2,3,4,5)
	ylocation=c(0.8,0.7,0.6,0.5)
	if(a==0){
		for(e in 1:4){
			abline(y0-x0*1/(alv[e]-1),1/(alv[e]-1),col=2,lwd=2)
			lab=paste("xi =1/",alv[e],sep="",collapse="")
			lscale=max(me)-y0
			a0=y0
			b0=1/(alv[e]-1)
			wy=y0+ylocation[e]*lscale
			wx=(wy-a0)/b0+x0
			text(wx,wy,lab)
		}
	}

	title(xlab="Threshold, u",ylab="Sample mean excess loss")
}



qqpareto=function(xv,pl=1,mle=0){
	n=length(xv)
	xv=xv[order(xv)]
	m1=mean(xv)
	m2=mean(xv^2)
	alpha=2*(m1*m1-m2)/(2*m1*m1-m2)
	lambda=m1*(alpha-1)
	if(mle==1){ pv=mle.pareto(xv); alpha=pv[1]; lambda=pv[2] }
	uv=((1:n)-0.5)/n
	xv3=qpareto(uv,lambda,alpha)
	if(pl == 1){
		plot(xv3,xv,xlab="",ylab="")
		title(xlab="Theoretical quantiles",ylab="Sample quantiles")
	}
	c(alpha,lambda)
}


mle.pareto=function(x,pr=1){
	pv0=c(4,1)
	pv=optim(pv0,ll.pareto,x=x,control=list(fnscale=-1))$par
	if(pr==1){
		cat("alpha  =",pv[1],"\n")
		cat("lambda =",pv[2],"\n")
		cat("log-likelihood =",ll.pareto(pv,x),"\n")
	}
	pv
}


ll.pareto=function(p,x){
	alpha=p[1]
	lambda=p[2]
	n=length(x)
	llv=log(dpareto(x,lambda,alpha))
	sum(llv)
}


ppareto=function(x,lambda,alpha){
	1-(lambda/(lambda+x))^alpha
}


qpareto=function(u,lambda,alpha){
	x=lambda/((1-u)^(1/alpha))-lambda
	x
}

rpareto=function(n,lambda,alpha){
	uv=runif(n)
	qpareto(uv,lambda,alpha)
}

dpareto=function(x,lambda,alpha){
	alpha*(lambda^alpha)/((lambda+x)^(alpha+1))
}



rgeometric=function(n,q){
	u=runif(n)
	x=floor(log(1-u)/log(q))+1
	x
}

qgeometric=function(u,q){
	x=floor(log(1-u)/log(q))+1
	x
}

pgeometric=function(x,q){
	p=1-q^x
	p
}

dgeometric=function(x,q){
	d=(1-q)*q^(x-1)
	d
}



qqgeometric=function(x){
        x=x[order(x)]	
        n=length(x)
        uv=((1:n)-0.5)/n
        q=(mean(x)-1)/mean(x)   
        xv=qgeometric(uv,q)	
        plot(xv,x,xlab="Theoretical quantiles",ylab="Sample quantiles",pch=20)
        title("QQ-Geometric Plot")
}



chi2test.geometric=function(x){
	n=length(x)
        q=(mean(x)-1)/mean(x) 
        m=5*floor(log(5/n)/log(q))
        ev=n*(1-q)*q^((1:m)-1)		
	C2threshold=5	
	bins=0
	lv0=1; lv1=1
	ec=0
	for(i in 1:m){
		ec=ec+ev[i]
		if(ec >= C2threshold){
			bins=bins+1
			lv1[bins]=i
			lv0[bins+1]=i+1
			ec=0
		}
	}
	lv1[bins]=m
	lv0=lv0[1:bins]

	nl=length(lv0)
	EV=(1:nl)*0
	OV=EV
	ov=ev*0
        for(e in 1:m){
                ov[e]=sum(x == e)
        }

	for(e in 1:nl){ EV[e]=sum(ev[lv0[e]:lv1[e]]); OV[e]=sum(ov[lv0[e]:lv1[e]]) }
	ch2=sum((OV-EV)^2/EV)
        df=length(OV)-2		
        pval=1-pchisq(ch2,df)
        cat("Time interval      \t Obs.\t Exp.\t Chi-sq\n")
        for(e in 1:(length(OV)-1)){
                cat(" t = ",lv0[e],"-",lv1[e],"\t\t",OV[e],"\t",round(EV[e],2),"\t",round((OV[e]-EV[e])^2/EV[e],3),"\n")
        }
                e=length(OV)
                cat(" t = ",lv0[e],"- oo","\t\t",OV[e],"\t",round(EV[e],2),"\t",round((OV[e]-EV[e])^2/EV[e],3),"\n")
        cat("Chi-squared statistic is ",ch2, "with ",df," degrees of freedom.\n")
        cat("The p-value for this is ",pval,".\n")
	cat("Minimum threshold for expected Bin frequency = ",C2threshold,".\n")
	


}


qqexp=function(xv,pl=1){
	n=length(xv)
	xv=xv[order(xv)]
	beta=1/mean(xv)
	uv=((1:n)-0.5)/n
	xv3=qexp(uv,beta)
	if(pl == 1){
		plot(xv3,xv,xlab="",ylab="")
		title(xlab="Theoretical quantiles",ylab="Sample quantiles")
		cat("beta = ",beta,"\n")
	}
	beta
}



ll.bt=function(p,x){
	x1=x[1,]
	x2=x[2,]
	mu1=p[1]
	mu2=p[2]
	mu=array(c(mu1,mu2),c(2,1))
	v11=p[3]
	v12=p[4]
	v21=v12
	v22=p[5]
	nu=p[6]
	detv=v11*v22-v12*v21
	va=array(c(v11,v21,v12,v22),c(2,2))
	vi=solve(va)
	k=-log(2*pi)-0.5*log(detv)
	llv=0
	llv=(x1-mu1)*vi[1,1]*(x1-mu1)+2*(x1-mu1)*vi[1,2]*(x2-mu2)+(x2-mu2)*vi[2,2]*(x2-mu2)
	llv=k-(nu+2)/2*log(1+1/nu*llv)
	ll=sum(llv)
	ll
}


mle.bt=function(x){
        options(warn=-1)   
	x1=x[1,]
	x2=x[2,]
	mu1=mean(x1)
	mu2=mean(x2)
	v11=var(x1)
	v22=var(x2)
	v12=mean((x1-mu1)*(x2-mu2))
	nu=10
	pv=c(mu1,mu2,v11,v12,v22,nu)
	pv1=optim(pv,ll.bt,x=x,control=list(fnscale=-1))
	cat(pv1$v,"\n")
	v1=-999999
	v0=pv1$v
	while(abs(v0-v1) > 0.01){
		v1=v0
		pv1=optim(pv1$par,ll.bt,x=x,control=list(fnscale=-1))
		cat(pv1$v,"\n")
		v0=pv1$v
	}
	options(warn=0)
	pv1$par
}

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
#'@importFrom graphics plot persp contour title
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
		setaxes(0,1,0,1)
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
		setaxes(0,1,0,1)
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
	
gcopula=
function(rho=0,l=50){
        uv1=((1:l)-0.5)/l
        uv2=uv1
        x1=qnorm(uv1)
        x2=x1
        da=array(0,c(l,l))
        for(i in 1:l){
                for(j in 1:l){
                        da[i,j]=dbnorm(c(x1[i],x2[j]),rho)
                }
        }
        ca=da*0
        for(i in 1:l){
                ca[i,]=cumsum(da[i,])
        }
        for(j in 1:l){
                ca[,j]=cumsum(ca[,j])
        }
        ca=ca/(l^2)
        ca2=array(0,c(l+1,l+1))
        ca2[2:(l+1),2:(l+1)]=ca
        u1=(0:l)/l
        u2=u1
        list(uv1=uv1,uv2=uv2,x1=x1,x2=x2,da=da,u1=u1,u2=u2,ca=ca2)
}


dbnorm=function(x,rho=0){
        x1=x[1]
        x2=x[2]
        detv=1-rho*rho
        density=-0.5*(x1*x1+x2*x2-2*rho*x1*x2)/detv
        density=1/(2*pi*sqrt(detv))*exp(density)/dnorm(x1)/dnorm(x2)
        density
}


tcopula=
function(df=4,rho=0,l=50){
        uv1=((1:l)-0.5)/l
        uv2=uv1
        x1=qt(uv1,df=df)
        x2=x1
        da=array(0,c(l,l))
        for(i in 1:l){
                for(j in 1:l){
                        da[i,j]=dbt(c(x1[i],x2[j]),df=df,rho=rho)
                }
        }
        ca=da*0
        for(i in 1:l){
                ca[i,]=cumsum(da[i,])
        }
        for(j in 1:l){
                ca[,j]=cumsum(ca[,j])
        }
        ca=ca/(l^2)
        ca2=array(0,c(l+1,l+1))
        ca2[2:(l+1),2:(l+1)]=ca
        u1=(0:l)/l
        u2=u1
        list(uv1=uv1,uv2=uv2,x1=x1,x2=x2,da=da,u1=u1,u2=u2,ca=ca2)
}

dbt=function(x,df=4,rho=0){
        x1=x[1]
        x2=x[2]
        detv=1-rho*rho
        g1=gamma(0.5*(df+2))
        g2=gamma(0.5*df)
        g3=(pi*df)
        density=1+(x1*x1+x2*x2-2*rho*x1*x2)/detv/df
        density=g1/g2/g3/sqrt(detv)*(density^(-(df+2)/2))/dt(x1,df=df)/dt(x2,df=df)
        density
}


plot.gauss.copula=function(rho=0,pl=1,zmax=10){
	res=gcopula(rho,l=100)
	if(pl==1){
		setaxes(0,1,0,1)
		vv=c((1:10)/10,1.5,2:10)
		title(xlab="u_1",ylab="u_2")
		contour(res$uv1,res$uv2,res$da,levels=vv,add=T)
	}
	else{
		nv=(0:33)*3+1
		persp(res$uv1[nv],res$uv2[nv],res$da[nv,nv],lwd=0.25,d=5,xlim=c(0,1),ylim=c(0,1),zlim=c(0,zmax),expand=0.75,xlab="u1",ylab="u2",zlab="copula density c(u1,u2)",ticktype="detailed",theta=300,phi=8)
	}
	a=paste("Gaussian Copula Density\nrho = ",round(rho,3),sep="",collapse="")
	title(a)
}





plot.t.copula=function(df=4,rho=0,pl=1,zmax=10){
	res=tcopula(df=df,rho=rho,l=100)
	if(pl==1){
		setaxes(0,1,0,1)
		vv=c((1:10)/10,1.5,2:10)
		title(xlab="u_1",ylab="u_2")
		contour(res$uv1,res$uv2,res$da,levels=vv,add=T)
	}
	else{
		nv=(0:33)*3+1
		persp(res$uv1[nv],res$uv2[nv],res$da[nv,nv],lwd=0.25,d=5,xlim=c(0,1),ylim=c(0,1),zlim=c(0,zmax),expand=0.75,xlab="u1",ylab="u2",zlab="copula density c(u1,u2)",ticktype="detailed",theta=300,phi=10)
	}
	a=paste("t Copula Density\ndf = ",round(df,2),"; rho = ",round(rho,3),sep="",collapse="")
	title(a)
}



ll.gumbel=function(theta,ua){
	u1=ua[1,]
	u2=ua[2,]
	x=(-log(u1))^theta+(-log(u2))^theta
	pd1=-theta/u1*((-log(u1))^(theta-1))
	pd2=-theta/u2*((-log(u2))^(theta-1))
	density=pd1*pd2*(1/theta^2)*((theta-1)*(x^(1/theta-2))+x^(2/theta-2))*exp(-x^(1/theta))
	log.density=log(density)
	sum(log.density)
}


mle.gumbel.copula=function(ua){
	theta1=optimize(ll.gumbel,c(1,50),ua=ua,maximum=TRUE)
	ll.max=theta1$obj
	theta1=theta1$max
	cat("Gumbel copula\n")
	cat("theta = ",theta1,"\n")
	cat("maximum likelihood = ",ll.max,"\n")
	theta1
}


ll.clayton=function(theta,ua){
	u1=ua[1,]
	u2=ua[2,]
	density=(theta+1)* ((u1*u2)^(-theta-1)) * ((u1^(-theta)+u2^(-theta)-1)^(-1/theta-2))
	log.density=log(density)
	sum(log.density)
}


mle.clayton.copula=function(ua){
	theta1=optimize(ll.clayton,c(0.01,50),ua=ua,maximum=TRUE)
	ll.max=theta1$obj
	theta1=theta1$max
	cat("Clayton copula\n")
	cat("theta = ",theta1,"\n")
	cat("maximum likelihood = ",ll.max,"\n")
	theta1
}


ll.gauss.copula=function(rho,ua){
	n=length(ua[1,])
	density=(1:n)*0
	for(i in 1:n){
		density[i]=dbnorm2(ua[,i],rho=rho)
	}
	log.density=log(density)
	sum(log.density)
}

dbnorm2=function(uv,rho=0){
	u1=uv[1]
	u2=uv[2]
	x1=qnorm(u1)
	x2=qnorm(u2)
	detv=1-rho*rho
	density=-0.5*(x1*x1+x2*x2-2*rho*x1*x2)/detv
	density=1/2/pi/sqrt(detv)*exp(density)/dnorm(x1)/dnorm(x2)
	density
}

mle.gauss.copula=function(ua){
	rho1=optimize(ll.gauss.copula,c(-0.999,0.999),ua=ua,maximum=TRUE)
	ll.max=rho1$obj
	rho1=rho1$max
	cat("Gaussian copula\n")
	cat("rho = ",rho1,"\n")
	cat("maximum likelihood = ",ll.max,"\n")
	rho1
}


ll.t.copula=function(pv,ua){
	nu=pv[1]
	rho=pv[2]
	n=length(ua[1,])
	density=(1:n)*0
	for(i in 1:n){
		density[i]=dbt2(ua[,i],df=nu,rho=rho)
	}
	log.density=log(density)
	sum(log.density)
}

dbt2=function(uv,df=4,rho=0){
	u1=uv[1]
	u2=uv[2]
	x1=qt(u1,df=df)
	x2=qt(u2,df=df)
	detv=1-rho*rho
	g1=0.5*df
	g2=1
	g3=(pi*df)
	density=1+(x1*x1+x2*x2-2*rho*x1*x2)/detv/df
	density=g1/g2/g3/sqrt(detv)*(density^(-(df+2)/2))/dt(x1,df=df)/dt(x2,df=df)
	density
}

mle.t.copula=function(ua){
	pv0=c(4,0.8)
	pv=optim(pv0,ll.t.copula,ua=ua,control=list(fnscale=-1))$par
	cat("t copula\n")
	cat("d.f. =",pv[1],"\n")
	cat("rho  =",pv[2],"\n")
	cat("log-likelihood =",ll.t.copula(pv,ua),"\n")
	pv
}


simulate.clayton.copula=function(theta,ns,m=1000){
	ns2=ns+sqrt(ns)*4
	u1=((1:m)-0.5)/m
	u2=u1
	v1=0
	v2=0	
	kk=0
	for(i in 1:m){
	if(i %% 10 == 0){ cat(i,"") }
	if(i %% 100 == 0){ cat("\n") }
		for(j in 1:m){
			zz=(theta+1)*((u1[i]*u2[j])^(-theta-1)) * ((u1[i]^(-theta)+u2[j]^(-theta)-1)^(-1/theta-2))
			n=rpois(1,zz*ns2/m/m)
			if(n > 0){
				for(k in 1:n){
					kk=kk+1
					uu1=(i-runif(1))/m
					uu2=(j-runif(1))/m
					v1[kk]=uu1
					v2[kk]=uu2
				}
			}
		}
	}
	cat("\n")
	va=t(array(c(v1,v2),c(length(v1),2)))
	if(kk > ns){
		uv=runif(kk)
		nv=order(uv)[1:ns]
		va=va[,nv]
	}
	va
}

simulate.gumbel.copula=function(theta,ns,m=1000){
	ns2=ns+sqrt(ns)*4
	u1=((1:m)-0.5)/m
	u2=u1
	v1=0
	v2=0	
	kk=0
	for(i in 1:m){
	if(i %% 10 == 0){ cat(i,"") }
	if(i %% 100 == 0){ cat("\n") }
		for(j in 1:m){
			vv1=u1[i]
			vv2=u2[j]
			x=(-log(vv1))^theta+(-log(vv2))^theta
			pd1=-theta/vv1*((-log(vv1))^(theta-1))
			pd2=-theta/vv2*((-log(vv2))^(theta-1))
			zz=pd1*pd2*(1/theta^2)*((theta-1)*(x^(1/theta-2))+x^(2/theta-2))*exp(-x^(1/theta))
			n=rpois(1,zz*ns2/m/m)
			if(n > 0){
				for(k in 1:n){
					kk=kk+1
					uu1=(i-runif(1))/m
					uu2=(j-runif(1))/m
					v1[kk]=uu1
					v2[kk]=uu2
				}
			}
		}
	}
	cat("\n")
	va=t(array(c(v1,v2),c(length(v1),2)))
	if(kk > ns){
		uv=runif(kk)
		nv=order(uv)[1:ns]
		va=va[,nv]
	}
	va
}


simulate.gauss.copula=function(rho,ns){
	z1=rnorm(ns)
	z2=rho*z1+sqrt(1-rho*rho)*rnorm(ns)
	u1=pnorm(z1)
	u2=pnorm(z2)
	va=t(array(c(u1,u2),c(ns,2)))
	va
}


simulate.t.copula=function(df=df,rho,ns){
	z1=rnorm(ns)
	z2=rho*z1+sqrt(1-rho*rho)*rnorm(ns)
	y=rchisq(ns,df=df)
	w1=z1/sqrt(y/df)
	w2=z2/sqrt(y/df)
	u1=pt(w1,df=df)
	u2=pt(w2,df=df)
	va=t(array(c(u1,u2),c(ns,2)))
	va
}



mle.nctB=function(x,nu,pr=1){
	pv0=c(mle.t(x,pr=0)[1:2],0)
	options(warn=-1)	
	pv=optim(pv0,ll.nctB,x=x,nu=nu,control=list(fnscale=-1))$par
	if(pr==1){
		cat("mu    =",pv[1],"\n")
		cat("sigma =",pv[2],"\n")
		cat("nu    =",nu,"\n")
		cat("ncp   =",pv[3],"\n")
		cat("log-likelihood =",ll.nctB(pv,x,nu),"\n")
	}
	options(warn=0)	
	c(pv[1],pv[2],nu,pv[3])
}

ll.nctB=function(p,x,nu){
	n=length(x)
	mu=p[1]
	sigma=p[2]
	ncp=p[3]
	y=(x-mu)/sigma
	llv=dt(y,df=nu,ncp=ncp,log=TRUE)-log(sigma)
	sum(llv)
}



chi2test.unif=function(u,constraints,pl=0){
	n=length(u)
	m=floor(n/10)
	vv=(0:m)/m
	Gv=vv*0
	for(e in 1:(m+1)){
		Gv[e]=sum(u < vv[e])
	}
	obs=diff(Gv)
	Ev=n/m
	chi2statistic=sum((obs-Ev)^2/Ev)
	pvalue=1-pchisq(chi2statistic,df=m-1-constraints)	
	if(pl == 1){ hist(u,breaks=vv) }
	cat("Test statistic, S = ",chi2statistic,"\n")
	cat("Number of observations = ",n,"\n")
	cat("Degrees of Freedom = ",m-1-constraints,"\n")
	cat("p-value = ",pvalue,"\n")
	pvalue
}


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

