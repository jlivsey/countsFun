#----------------------------------------------------------------------------------#
# Adding here some functions from the now obsolete MixtureInf package.
# https://cran.r-project.org/src/contrib/Archive/MixtureInf/
#
# These may be used in the initial estimation of Mixed Poisson
#----------------------------------------------------------------------------------#
# Combine round and signif functions.
rousignif = function(x){
  #x:	a single numeric or a numeric vector or matrix.
    a=NULL
    if (is.matrix(x))
    {
      a=c(nrow(x),ncol(x))
      b1=rownames(x)
      b2=colnames(x)
      x=as.vector(x)
    }
    x[abs(x)>=0.001]=round(x[abs(x)>0.001],3)
    x[abs(x)<=0.001]=signif(x[abs(x)<=0.001],3)
    if (is.null(a)==F)
    {
      x=matrix(x,a)
      rownames(x)=b1
      colnames(x)=b2
    }
    x
  }

# Compute the PMLE of parameters under a Poisson mixture
pmle.pois = function(x,m0=1,lambda=0,inival=NULL,len=10,niter=50,tol=1e-6,rformat=FALSE){
    #x: 		data, can be either a vector or a matrix with the 1st column being the observed values
    #        		and the 2nd column being the corresponding frequencies.
    #m0: 		order of finite mixture model.
    #inival:	initial values for the EM-algorithm
    #lambda:	size of penalty function of mixing proportions.
    #len: 	number of initial values chosen for the EM-algorithm.
    #niter:     least number of iterations for all initial values in the EM-algorithm.
    #tol: 	tolerance value for the convergence of the EM-algorithm.
    #rformat	format for output, rformat=T means the format of output is determined by R software.
    #		rformat=F means the format of output is determined by our default setting. When the output is
    #		larger than 0.001, it is determined by round(output,3); When the output is less than 0.001,
    #		it is determined by signif(output,3).

    if (is.data.frame(x))
    {
      if (ncol(x)==2)
        x=as.matrix(x)
      if (ncol(x)==1 | ncol(x)>2)
        x=x[,1]
    }
    if(is.vector(x))
    {
      y=as.matrix(table(x))
      count=as.numeric(rownames(y))
      freq=y[,1]
      x=cbind(count,freq)
    }

    out=phi0.pois(x,m0,lambda,inival,len,niter,tol)
    alpha=out$alpha
    theta=out$theta
    loglik=out$loglik
    ploglik=out$ploglik

    if (rformat==F)
    {
      alpha=rousignif(alpha)
      theta=rousignif(theta)
      loglik=rousignif(loglik)
      ploglik=rousignif(ploglik)
    }

    if (lambda==0 | m0==1)
      list('MLE of mixing proportions:'=alpha,
           'MLE of component parameters:'=theta,
           'log-likelihood:'=loglik)
    else
      list('PMLE of mixing proportions:'=alpha,
           'PMLE of component parameters:'=theta,
           'log-likelihood:'=loglik,
           'Penalized log-likelihood:'=ploglik)
  }


# Compute the PMLE of mixing distribution under null hypothesis for Poisson mixture.
phi0.pois = function(x,m0,lambda,inival,len,niter,tol){
    #x:      data; can be either a vector or a matrix with the 1st column being the observed values
    #        		and the 2nd column being the corresponding frequencies.
    #m0:     order under null hypothesis.
    #lambda:	size of penalty function of mixing proportions.
    #inival: initial values for the EM-algorithm.
    #len: 	number of initial values chosen for the EM-algorithm.
    #niter:     least number of iterations for all initial values in the EM-algorithm.
    #tol: 	tolerance value for the convergence of the EM-algorithm.

    count=as.numeric(x[,1])
    freq=as.numeric(x[,2])
    if(m0>1)
    {
      if (is.vector(inival))
        inival=t(inival)
      if (is.null(inival)==F)
        len=nrow(inival)
      output=c()
      for(i in 1:len)
      {
        if (is.null(inival))
        {
          alpha=runif(m0,0,1)
          alpha=alpha/sum(alpha)
          theta=sort(runif(m0,0,max(count)))
        }
        else
        {
          alpha=inival[i,1:m0]
          alpha=alpha/sum(alpha)
          theta=sort(inival[i,(m0+1):(2*m0)])
        }
        for (j in 1:niter)###run niter EM-iterations first
        {
          pdf.component=t(t(apply(as.matrix(theta,ncol=1),1,dpois,x=count))*alpha)+1e-100/m0
          pdf=apply(pdf.component,1,sum)
          w=pdf.component/pdf
          alpha=(apply(freq*w,2,sum)+lambda)/(sum(freq)+m0*lambda)
          theta=apply(freq*w*count,2,sum)/apply(freq*w,2,sum)
        }
        pdf.component=t(t(apply(as.matrix(theta,ncol=1),1,dpois,x=count))*alpha)+1e-100/m0
        pdf=apply(pdf.component,1,sum)
        pln=sum(freq*log(pdf))+lambda*sum(log(alpha))
        output=rbind(output,c(alpha,theta,pln))
      }
      index=which.max(output[,(2*m0+1)])
      alpha=output[index,1:m0]
      theta=output[index,(m0+1):(2*m0)]
      pln0=output[index,(2*m0+1)]
      err=1
      t=0
      pdf.component=t(t(apply(as.matrix(theta,ncol=1),1,dpois,x=count))*alpha)+1e-100/m0
      pdf=apply(pdf.component,1,sum)
      while(err>tol & t<2000)###EM-iteration with the initial value with the largest penalized log-likelihood
      {
        w=pdf.component/pdf
        alpha=(apply(freq*w,2,sum)+lambda)/(sum(freq)+m0*lambda)
        theta=apply(freq*w*count,2,sum)/apply(freq*w,2,sum)
        pdf.component=t(t(apply(as.matrix(theta,ncol=1),1,dpois,x=count))*alpha)+1e-100/m0
        pdf=apply(pdf.component,1,sum)
        pln1=sum(freq*log(pdf))+lambda*sum(log(alpha))
        err=abs(pln1-pln0)
        pln0=pln1
        t=t+1
      }
      ln=pln1-lambda*sum(log(alpha))
      index=sort(theta,index.return=TRUE)$ix
      alpha0=alpha[index]
      theta0=theta[index]

      list("alpha"=alpha0,"theta"=theta0,"loglik"=ln,"ploglik"=pln1)
    }
    else
    {
      theta0=sum(freq*count)/sum(freq)
      list("alpha"=1,"theta"=theta0,"loglik"=sum(freq*log(dpois(count,theta0))),
           "ploglik"=sum(freq*log(dpois(count,theta0))))
    }
  }
