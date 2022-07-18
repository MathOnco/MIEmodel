## ODE Model equation
chrmod <- function(time,state,parms){
  with(as.list(parms),{
    ds <- state%*%A
    return(list(ds))
  })
}

run_ode <- function(beta,mu, maxchrom, founder, minchrom=1,times=0:1000,
                    B=function(beta,i) return(beta),
                    D=function(mu,i) return(mu)){
  Ndim <- length(maxchrom)
  Nstates <- prod(maxchrom-minchrom+1)
  state <- rep(0,Nstates)
  state[state_index(founder,maxchrom,minchrom)] <- 1
  names(state) <- paste0("(",sapply(1:Nstates, function(i){
    paste(reverse_state_index(i,maxchrom,minchrom,Ndim),
          collapse=",")
  }),")")
  A <- get_A(maxchrom,B,D,beta,mu,minchrom)
  parms <- list(A=A)
  out <- ode(y=state,func=chrmod,parms=parms,times = times)
  out <- out[!is.na(rowSums(out)),]
  data.frame(out,check.names = F)
}

## returns true if population grew, false otherwise
check_growth <- function(beta,mu, maxchrom, founder, minchrom=1, times=0:1000,
                 B=function(beta,i) return(beta),
                 D=function(mu,i) return(mu)){
  
  out <- run_ode(beta,mu, maxchrom, founder, minchrom, times, B,D)
  sum(out[1,-1])<sum(out[nrow(out),-1])
}

## Define average number of cells with j chromosomes resulting from one cell with i chromosomes in the previous generation
pij<-function(i, j, beta){
  qij <- 0
  if(abs(i-j)>i){ ## not enough copies for i->j
    return(qij)
  }
  # code fails for j = 0, but we can get this result by noting i->0 always
  ## accompanies i->2i
  if(j==0) j <- 2*i
  s <- seq(abs(i-j), i, by=2 )
  for(z in s){
    qij <- qij + choose(i,z) * beta^z*(1-beta)^(i-z) * 0.5^z*choose(z, (z+i-j)/2)
  }
  ## if there is a mis-segregation: ploidy conservation implies that both daughter cells will emerge simultaneously
  #if(i!=j) qij <- 2*qij
  
  return(qij)
}

# maps copy number state to a unique integer
state_index <- function(state,maxchrom,minchrom=1){
  if(prod(minchrom<= state & state<=maxchrom )==0) stop("Error: invalid state")
  ndim <- length(state)
  mc <- maxchrom-minchrom+1
  if(length(mc)==1 & length(state)>1) mc <- rep(mc,length(state))
  sx <- state-minchrom+1
  cp <- c(1,cumprod(mc))
  sum(sapply(1:length(state), function(i) (sx[i]-1)*cp[i]))+1
}

# reverses mapping from integer to state
reverse_state_index <- function(index,maxchrom,minchrom=1,ndim){
  ## reset maxchrom to behave as if minchrom was 1,1,...
  mc <- maxchrom-minchrom+1
  if(length(mc)==1 & ndim>1) mc <- rep(mc,ndim)
  ## how many elements does each part of the state represent?
  ## works as prod(numeric(0)) evaluates to 1:
  Nsites <- cumprod(mc)
  cp <- c(1,cumprod(mc)[-length(mc)])
  state <- c()
  for(j in ndim:1){
    ni <- floor((index-1)/cp[j])
    state <- c(ni+1,state)
    index <- index-cp[j]*ni
  }
  state + minchrom-1
}

#reverse_state_index(state_index(c(2,2,2),maxchrom=8,minchrom=3),maxchrom=8, minchrom=2,ndim=3)

get_A <- function(maxchrom,B=function(beta,i) return(beta),
                  D=function(mu,i) return(mu), beta,mu,minchrom=1){
  
  Ndim <- length(maxchrom)
  Nstates <- prod(maxchrom-minchrom+1) 
  A<- matrix(0,nrow=Nstates,ncol=Nstates)
  for(i in 1:nrow(A)){
    state_i <- reverse_state_index(i,maxchrom,minchrom,ndim=Ndim)
    beta_i <- B(beta,state_i)
    mu_i <- D(mu,state_i)
    #beta_i <- B(beta,i)
    #mu_i <- D(mu,i)  
    for(j in 1:nrow(A)){
      state_j <- reverse_state_index(j,maxchrom,minchrom,ndim=Ndim)
      
      ##individual transition probabilities:
      qij <- sapply(1:Ndim, function(k) pij(state_i[k], state_j[k], beta_i))
      
      ## joint probability (i1,i2,...)->(j1,j2,...)
      qij <- prod(qij)
      A[i,j] <- 2*qij
      
      # ## case when there is no mis-segregation:
      if(i==j) A[i,j] <- (2*qij-1)-mu_i
      
    }
  }
  A
}

## determine whether given parameters lie on the critical curve (occurs when function returns 0)
c_crit <- function(mu,beta,maxchrom, minchrom=1,
                   B=function(beta,i) return(beta),
                   D=function(mu,i) return(mu)){
  

  ## generate transition matrix
  #A <- get_A(maxchrom,beta,mu,minchrom,B,D)
  A <- get_A(maxchrom,B,D,beta,mu,minchrom)
  ## find eigenvalues
  z <- tryCatch(eigen(A)$value, error = function(e) return("error"))
  if(z[1]=="error") return(10^5)
  res <- 10^5
  if(sum(Im(z)==0)>0){
    res <- abs(max(Re(z[Im(z)==0]))) ##distance of the real part from zero
  }
  return(res)
}

## find the value of mu lying closest to the critical curve for a specified beta
find_c_crit <- function(beta,maxchrom,minchrom=1, mu_range = c(0,1),
                        B=function(beta,i) return(beta),
                        D=function(mu,i) return(mu)){

  out <- optimise(c_crit,interval=mu_range,tol=0.000000000000000001,beta=beta,B=B,maxchrom=maxchrom,minchrom=minchrom,D=D)
  res <- data.frame(v=out$objective,mu=out$minimum,beta=beta,maxchrom=paste(maxchrom,collapse=","))
  return(res)
}

## following function works for multiple chromosomes if maxchrom is a vector. 
## However, still limited because for n chromosomes we will have to generate 
## a matrix with size maxchrom^n
#v            mu          beta   maxchrom
#4.072862e-09 1.205202    2      110
fast_c_crit <- function(beta,maxchrom,minchrom=1,
                        B=function(beta,i) return(beta),
                        D=function(mu,i) return(mu),df=FALSE){
  ## generate critical curve setting mu=0
  A <- get_A(beta=beta,maxchrom=maxchrom,minchrom=minchrom,mu=0,B=B,D=D)
  if(!df) return(max(Re(eigen(A)$values)))
  return(return(data.frame(v=0,mu=Re(max(Re(eigen(A)$values))),beta=beta,maxchrom=paste(maxchrom,collapse=","))))
  ## find steady state distribution - note this is invariant to mu
  #ss <- eigen(t(A))
  #ss <- ss$vectors[,Re(ss$values)==max(Re(ss$values))]
  #ss <- ss/sum(ss)
  ## calculate the death flux due to mis-segregations only
  ## (the reason we set mu=0 was to make this easier)
  #f_flux <- 1-rowSums(A)
  #dx <- sum(ss*f_flux)
  ## return value of mu on the critical curve
  #if(!df) return(Re(1-dx))
  #return(data.frame(v=0,mu=Re(1-dx),beta=beta,maxchrom=paste(maxchrom,collapse=",")))
}

## following function extrapolates from one chromosome to multiple chromosomes.
## Could be further adapted to handle chromosomes with different bounds. 
very_fast_c_crit <- function(beta,maxchrom,n,minchrom=1,
                             B=function(beta,i) return(beta),
                             D=function(mu,i) return(mu)){
  ## generate critical curve setting mu=0
  A <- get_A(beta=beta,maxchrom=maxchrom,minchrom=minchrom,mu=0,B=B)
  ## recover the matrix qij from A (could have just built this matrix directly I guess)
  q <- do.call(rbind,lapply(1:nrow(A), function(i){
    sapply(1:ncol(A), function(j){
      qij <- A[i,j]/2
      if(i==j) qij <- (A[i,j]+1)/2
      return(qij)
    })
  }))
  
  ## these are the probabilities that each cell with cn=i mis-segregates
  ## to a dead state:
  pxn <- 1-rowSums(q)
  
  ## calculate the steady state distribution
  ss <- eigen(t(A))
  ss <- ss$vectors[,Re(ss$values)==max(Re(ss$values))]
  ss <- ss/sum(ss)
  
  ## (1-sum(ss*pxn))^n gives the probability that a cell at steady state with n chromosomes
  ## does not die from mis-segregation during division
  ## We can get the flux from the ODE equation, i.e. where we had:
  ## A[i,j] <- (2*qij-1)-mu_i
  ## at the critical curve, we want to set the net rate of change of the 
  ## alive population = 0, effectively:
  ## A[alive,alive]=0
  ## so 0 = 2*qij-1 - mu_ccrit
  mu_ccrit <- Re(2*(1-sum(ss*pxn))^n-1)
  return(mu_ccrit)
}

pop_ss_B <- function(beta,mu,minchrom=1,maxchrom,B=function(beta,i) return(beta),
                     D=function(mu,i) return(mu)){
  A <- get_A(maxchrom=maxchrom,minchrom=minchrom,beta=beta,mu=mu,B=B,D=D)
  ss <- eigen(t(A))
  ss <- ss$vectors[,Re(ss$values)==max(Re(ss$values))]
  ss <- ss/sum(ss)
  p_mis <- sapply(minchrom:maxchrom, function(i) B(beta,i))
  sum(ss*p_mis)
}


pop_ss_D <- function(beta,mu,maxchrom,minchrom=1,B=function(beta,i) return(beta),
                     D=function(mu,i) return(mu)){
  A <- get_A(maxchrom=maxchrom,minchrom=minchrom,beta=beta,mu=mu,B=B,D=D)
  ss <- eigen(t(A))
  ss <- ss$vectors[,Re(ss$values)==max(Re(ss$values))]
  p_D <- sapply(minchrom:maxchrom, function(i) D(mu,i))
  #ss <- ss[p_D<1]
  #p_D <- p_D[p_D<1]
  ss <- (ss/sum(ss))
  sum(ss*p_D)
}



boxplot2D <- function(df, critcurve, pdfout){
  # ---
  # author: "Thomas Veith"
  # date: "01/02/2021"
  # ---
  
  xlab = expression(paste("Turnover rate (",mu,"/",lambda, ")"))
  ylab = expression(paste("Mis-segregation Rate (", beta, ")"))
  col=c("cyan","orange","purple")
  names(col)=c("Breast","Lung","Skin")
  
  # pdf(pdfout, width = 5,height = 3.3)
  p <- ggplot(df, aes(fill = group, color = group)) +
    geom_line(
      data = critcurve,
      aes(x = mu, y = beta),
      size = 1,
      color = "black",
      linetype = "solid"
    )  +
    
    ## 2D box defined by the Q1 & Q3 values in each dimension, with outline
    geom_rect(aes(
      xmin = x.lower,
      xmax = x.upper,
      ymin = y.lower,
      ymax = y.upper
    ),
    alpha = 0.3) +
    geom_rect(
      aes(
        xmin = x.lower,
        xmax = x.upper,
        ymin = y.lower,
        ymax = y.upper
      ),
      color = "black",
      fill = NA
    ) +
    
    ## whiskers for x-axis dimension with ends
    geom_segment(aes(
      x = x.min,
      y = y.middle,
      xend = x.max,
      yend = y.middle
    )) + ## whiskers
    geom_segment(aes(
      x = x.min,
      y = y.lower,
      xend = x.min,
      yend = y.upper
    )) + ## lower end
    geom_segment(aes(
      x = x.max,
      y = y.lower,
      xend = x.max,
      yend = y.upper
    )) + ## upper end
    
    ## whiskers for y-axis dimension with ends
    geom_segment(aes(
      x = x.middle,
      y = y.min,
      xend = x.middle,
      yend = y.max
    )) + ## whiskers
    geom_segment(aes(
      x = x.lower,
      y = y.min,
      xend = x.upper,
      yend = y.min
    )) + ## lower end
    geom_segment(aes(
      x = x.lower,
      y = y.max,
      xend = x.upper,
      yend = y.max
    )) + ## upper end
    
    ylab(ylab) + xlab(xlab)+
    scale_y_continuous(trans = "log", breaks=c(0, 0.002, 0.05, 1))  +
    theme_classic() +
    scale_fill_manual(values=c(col[df$group],"black")) +
    scale_color_manual(values=c(col[df$group],"black")) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(colour = "black", size = 4)
    ) 
  
  ggsave(pdfout,width = 5,height = 3.3)
}



