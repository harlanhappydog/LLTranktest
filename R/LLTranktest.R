#' @title LLT rank-based non-parameteric test
#'
#' @description LLT rank-based non-parameteric test
#'
#' @param formula, id=NULL, data, B
#'
#' @return data.frame
#'
#' @examples LLTranktest(UNL_mm ~ month + treat_group, id=toenail$ID, data=toenail, B=21)
#'
#' @export LLTranktest

####################################################
LLTranktest<-function(formula, id=NULL, data, B=200, Q=10){

  require(SphericalCubature)
  ##### The setup
  cl <- match.call()

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1,m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- quote(stats::model.frame)
  mf <- eval.parent(mf)

  mt <- attr(mf, "terms")
  X <- model.matrix(mt, mf, contrasts)[,-1]
  Y <- model.response(mf, 'numeric')

  if(is.null(id)){id<-paste(1:length(Y))}

  X<-as.matrix(cbind(X))
  p<-dim(X)[2]
  N<-length(unique(id))
  NM<-length(Y)
  s<-rank(Y)

  if(p==1){warning("Test only available for more than one covariate")
    break}

########################


  ## Using (p-1) dimensional polar coords
sn_pdim<-function(mythetavec, psi=rep(1,N), h_mult=1){

    h= h_mult*(N^(-1/3))
    reppsi<-c(unlist(apply(cbind(1:N),1,function(x) rep(psi[x],sum(id==unique(id)[x])))))
    betavec<-cbind(rev(c(polar2rect(1, mythetavec))))

    okaym<-function(m){ sum(reppsi*reppsi[m]*(s>s[m])*pnorm((X%*%betavec-c(X[m,]%*%betavec))/h))}
    bb<-sum(-apply(cbind(1:NM),1,okaym))/((N)*(N-1))
    return(bb)
  }
  ##########################
  
 generateoriginal<-function(my_H){
  
  sn_pdim_H<-function(x){sn_pdim(x, psi=rep(1,N), h_mult=my_H)}
  
  init_points<-matrix(runif(6*(p-1), -3.142,  3.142),6,)
  trial_outs<-init_points*0
  trials<-rep(0,6)
  for(jj in 1:dim(init_points)[1]){
  	
  	tryouts<-nlm(sn_pdim_H, init_points[jj,], iterlim=3)
  	trial_outs[jj,]<-tryouts$estimate
  	trials[jj]<-tryouts$minimum
  	
  }
  theta_hat = nlm(sn_pdim_H, trial_outs[which.min(trials),])$estimate
  mybetavec<-rev(polar2rect(1,theta_hat))
  yfit = c(X%*%mybetavec)
  return(list(sd=sd(yfit), point_est=mybetavec))}
  
  ##########################

myH<-generateoriginal(1)$sd
print("Part 1: establishing the bandwidth and point estimate.")
for(q in 1:Q){
part1_out<-generateoriginal(myH)
myH_new<-part1_out$sd
myH<-myH_new
print(paste("sigma[",q,"]=", round(myH,3), sep=""))}
my_betavec<-part1_out$point_est
print(paste("point estimate ="))
print(round(my_betavec,4))

print("Part 2: Calculating p-values.")
  mybetavec_b<-matrix(0,B,p)
  pb = txtProgressBar(min = 0, max = B, initial = 0)
  for(bb in 1:B){
    setTxtProgressBar(pb,bb)

    tryCatch({
      psib=rexp(N)
      snB_pdim<-function(x){sn_pdim(mythetavec=x, psi=psib, h_mult= myH)}

        init_points<-matrix(runif(6*(p-1), -3.142,  3.142),6,)
  trial_outs<-init_points*0
  trials<-rep(0,6)
  for(jj in 1:dim(init_points)[1]){
  	
  	tryouts<-nlm(snB_pdim, init_points[jj,], iterlim=3)
  	trial_outs[jj,]<-tryouts$estimate
  	trials[jj]<-tryouts$minimum
  	
  }
      theta_hat_b = nlm(snB_pdim, trial_outs[which.min(trials),])$estimate
      mybetavec_b[bb,]<-rev(polar2rect(1, theta_hat_b))


    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  close(pb)

  pval1 <-rep(0,p)
  for(pp in 1:dim(X)[2]){
    pval1[pp]<-	2*min(c(   (1+ sum(mybetavec_b[,pp]>0))/(1+B) ,  (1+ sum(mybetavec_b[,pp]<=0))/(1+B)  ))
  }

  CI<-apply(mybetavec_b,2,function(x) quantile(x,c(0.025,0.975)))
  result<-data.frame(scaled_estimate=c(my_betavec), lower95CI=CI[1,],upper95CI=CI[2,], pval= pval1)
  rownames(result)<-colnames(X)

  print(paste("Number of Observations:",NM))
  print(paste("Number of Groups:",N))

  return(result)}

