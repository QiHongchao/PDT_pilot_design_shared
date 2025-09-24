##design matrix
sw.design.mat <- function(I,J,H=NULL) {
  ## Creates the design matrix for a stepped wedge design
  ## Checks to see that data are consistent with SW design
  if(sum(sapply(list(I,J,H), is.null)) != 1) {
    warning("exactly one of 'I', 'J' and 'H' must be NULL")
  }
  if (is.null(I)) {
    I <- H*J
  }
  if (is.null(J)) {
    J <- I/H
  }
  if (is.null(H)) {
    H <- I/J
  }
  
  # Stepped wedge design matrix
  X <- matrix(0,I,(J+1))            
  for (i in 2:(J+1)) {
    X[1:((i-1)*H),i] <- 1
  }
  row.names(X) <- sample(1:I,I)
  colnames(X) <- c("Baseline",paste0("Time ",1:J))
  return(X)
}
##simulate a trial
make.swt <- function(I=NULL,J=NULL,H=NULL,K,design="cross-sec",mu=NULL,b.trt,
                     b.time=NULL,sigma.y=NULL,sigma.e=NULL,rho,sigma.a=NULL,rho.ind=NULL,
                     sigma.v=NULL,X=NULL,family="binomial",natural.scale=TRUE)  {
  ## Normal outcome
  ## inputs:
  # I = number of clusters
  # J = number of randomisation time points (excluding baseline)
  # K = average sample size in each cluster
  # H = number of units randomised at each time point
  # design = type of design. Can be "cross-sec" (default) or "cohort" (repeated measurements)
  # mu = baseline mean
  # b.trt = intervention effect
  # b.time = (optional) time effect (on linear time trend scale!)
  # sigma.y = *total* individual variability
  # sigma.e = *residual* individual variability
  # rho = ICC at the cluster level
  # sigma.a = the sd of the cluster random effect (default at NULL)
  # rho.ind = ICC at the individual level (for cohort design, default at NULL)
  # sigma.v = the sd of the treatment random effect (default at NULL)
  # X = SW design matrix (default at NULL and will be computed automatically)
  # family = type of outcome to be simulated (options are "gaussian", "binomial" and "poisson")
  # natural.scale = whether the input values for the trt effect etc are passed on the natural scale or not
  #
  # GB + Rosalind Leach (November 2015)
  
  # CHECK ON VALID FAMILY TYPE 
  flag=1
  if(family=="gaussian" | family=="binomial" | family=="poisson"){flag=0}
  if(flag==1){ stop("Available options for the argument 'family' are: gaussian, binomial or poisson")}
  
  # CREATE DESIGN MATRIX FROM sw.design.mat IF X HAS NOT BEEN PROVIDED
  if(is.null(X)) {X <- sw.design.mat(I=I,J=J,H=NULL)} else {row.names(X) <- sample(1:I,I)}
  
  
  if(family=="binomial"){
    # SETS DEFAULT FOR mu AND b.time IF THEY HAVE NOT BEEN PORVIDED
    if(is.null(mu)) {
      if (natural.scale==TRUE) {
        mu=.5		}		      # Default baseline probability on [0-1]
      else if (natural.scale==FALSE) {
        mu=log(.5/.5)	}		# Default baseline probability on logit scale
    } 
    
    # SIMULATES VARIABLES FROM THE INPUTS GIVEN
    treatment <- rep(t(X),each=K,1)                 # Treatment allocation at each time point
    time <- rep(seq(0,J),each=K,I)                  # K measurements during the trial
    person <- seq(1,K)                              # individuals IDs
    cluster.order <- rep(sample(1:I))               # order at which clusters switch
    cluster <- rep(cluster.order, each=(J+1)*K)   	# cluster IDs (randomised)
    id <- NULL
    
    # If natural.scale = T, then assume that the user is giving values for p0 (mu) and OR directly (b.trt)
    if (natural.scale==TRUE) {
      OR <- b.trt			              # for simplicity defines the OR
      b.trt <- log(b.trt)			      # logOR to be used in the linear predictor
      p0 <- mu				              # for simplicity defines the baseline probability
      mu <- log(mu/(1-mu))          # Converts baseline probability to log odds to be used in the linear predictor
    }  
    # But if natural.scale = F, then the user has passed data on the logit scale for p0 (mu) and on the logOR (b.trt)
    if (natural.scale==FALSE) {
      OR <- exp(b.trt)			        # Defines the OR on the natural scale
      p0 <- exp(mu)/(1+exp(mu))		  # Rescales the baseline to the [0-1] range
    }
    # Now defines p1 and sigma.e
    p1 <- OR*(p0/(1-p0))/(1+OR*(p0/(1-p0)))	      # Estimates the outcome probability for intervention
    sigma.e <- sqrt(((p0*(1-p0))+(p1*(1-p1)))/2)	# Pooled estimate of the within cluster sd
    if(is.null(b.time)) {                         # Time effect (if not specified, set to)
      b.time <- .5*b.trt  }                       #   some pre-defined proportion of the treatment effect
    
    # CALCULATE THE HYPERPARAMETERS
    mu.a <- 0		                              # Within cluster (random intercept) mean 
    
    # CROSS-SECTIONAL DESIGN -------------------------------------------------
    if(design=="cross-sec") {
      # CHECK THE RIGHT COMBINATION OF sigma.e, sigma.y, sigma.a, rho HAS BEEN GIVEN
      if((sum(sapply(list(sigma.y,sigma.e),is.null)) !=1) & (sum(sapply(list(sigma.a,rho),is.null)) !=1)) { 	# Need to pass 1 of 'sigma.y' or 'sigma.e' and 1 of 'sigma.a' or 'rho
        stop("Please provide either 'sigma.y' (*total* individual variability) or 'sigma.e' (*residual* variability)
             and either 'sigma.a' (*within cluster* variability) or 'rho' (ICC)")}
      
      # IF STATEMENTS TO CALCULATE ALL RELEVANT PARAMETERS FROM THE INPUTS GIVEN
      if(is.null(sigma.e) & is.null(sigma.a)) {
        sigma.a <- sqrt(rho*sigma.y^2)
        sigma.e <- sqrt(sigma.y^2-sigma.a^2)
      }
      if(is.null(sigma.e) & is.null(rho)) {
        sigma.e <- sqrt(sigma.y^2-sigma.a^2)
        rho <- sigma.a^2/sigma.y^2
      }
      ##the situation we have!!!!!!!
      if(is.null(sigma.y) & is.null(sigma.a)) {
        sigma.a <- sqrt(sigma.e^2*rho/(1-rho))
        sigma.y <- sqrt(sigma.a^2+sigma.e^2)
      }
      if(is.null(sigma.y) & is.null(rho)) {
        sigma.y <- sqrt(sigma.e^2+sigma.a^2)
        rho <- sigma.a^2/sigma.y^2
      }
      
      # SIMULATE REMAINING VARIABLES
      a <- rnorm(I, mu.a, sigma.a)  		            # cluster-level intercept
      
      # CALCULATE LINEAR PREDICTOR
      linpred <- mu + a[cluster] + b.trt*treatment + b.time*time		# Calculates the linear predictor
    } 
    
    # COHORT DESIGN -------------------------------------------------
    # if (design=="cohort") {  
    #   # CHECK THE RIGHT COMBINATION OF sigma.a, sigma.v, rho, rho.ind HAS BEEN GIVEN
    #   if((sum(sapply(list(rho, sigma.a),is.null)) !=1) || (sum(sapply(list(sigma.a,rho),is.null)) !=1) ||
    #      (sum(sapply(list(sigma.v,rho.ind),is.null)) !=1)) {
    #     stop("Please provide either 'sigma.a' (*within cluster* variability) or 'rho' (cluster-level ICC)
    #          AND either 'sigma.v' (*within cluster-individual* variability) or 'rho.ind (individual-level ICC)") }
    #   
    #   # IF STATEMENTS TO CALCULATE ALL RELEVANT PARAMETERS FROM THE INPUTS GIVEN
    #   if(is.null(sigma.v) & is.null(sigma.a)) {
    #     r1 <- rho.ind/(1-rho.ind)
    #     r2 <- rho/(1-rho)
    #     sigma.v <- sqrt((r1*(sigma.e^2/(1-rho)))/(1-r1*r2))
    #     sigma.a <- sqrt(r2*(sigma.e^2+sigma.v^2))
    #   }
    #   
    #   if(is.null(rho.ind) & is.null(sigma.a)) {
    #     r <- rho/(1-rho)
    #     sigma.a <- sqrt(r*(sigma.v^2+sigma.e^2))
    #     #rho.ind <-
    #   }
    #   
    #   if(is.null(rho.ind) & is.null(rho)) {
    #     rho <- sigma.a^2/(sigma.a^2+sigma.e^2+sigma.v^2)
    #     rho.ind <- sigma.v^2/(sigma.a^2+sigma.e^2+sigma.v^2)
    #   }
    #   
    #   if(is.null(rho) & is.null(sigma.v)) {
    #     r <- rho.ind/(1-rho.ind)
    #     sigma.v <- sqrt(r*(sigma.a^2+sigma.e^2))
    #     rho <- sigma.a^2/(sigma.a^2+sigma.e^2+sigma.v^2)
    #   }
    #   
    #   # SIMULATE REMAINING VARIABLES
    #   start <- seq(1,I*K,by=K)
    #   id <- numeric                                 # Individual ID (for repeated measurements)
    #   id <- rep(seq(start[cluster.order[1]],start[cluster.order[1]]+K-1),(J+1))
    #   for (i in 2:I) {
    #     id <- c(id,rep(seq(start[cluster.order[i]],start[cluster.order[i]]+K-1),(J+1)))
    #   }
    #   mu.v <- 0                                     # within individuals (random intercept) mean
    #   v <- rnorm(K*I,mu.v,sigma.v)                  # individual-level intercept
    #   ###      sigma.a <- sqrt(sigma.y^2*rho)                # Within cluster (random intercept) sd 
    #   a <- rnorm (I, mu.a, sigma.a)                 # cluster-level intercept 
    #   
    #   # CALCULATE LINEAR PREDICTOR ON THE LOGIT SCALE!
    #   linpred <- mu + a[cluster] + v[id] + b.trt*treatment + b.time*time
    # } 
    
    # SIMULATES DATA 
    n <- I*(J+1)*K
    # Converts log odds back to probability
    p = (exp(linpred))/(1+exp(linpred))
    y <- rbinom(n, 1, p)
  }
  
  # RETURNS DATA FRAME
  data <- data.frame(y, person, time, cluster, treatment, linpred, b.trt)
  if (design=="cohort") {data <- data.frame(data,id)}
  return(data)
}

##calculate the power with a new function
sim.power.ni.simplified <- function (I,J,H=NULL,K,design="cross-sec",mu=0,b.trt,b.time=NULL,
                           sigma.y=NULL,sigma.e=NULL,rho=NULL,sigma.a=NULL,
                           rho.ind=NULL,sigma.v=NULL,n.sims=1000,formula=NULL,
                           family="binomal",natural.scale=TRUE,sig.level=0.05,
                           method="lme",ni.margin,...) {
  
  exArgs <- list(...)
  
  requireNamespace("lme4",quietly=TRUE)
  requireNamespace("foreach",quietly=TRUE)
  requireNamespace("doParallel",quietly=TRUE)
  requireNamespace("parallel",quietly=TRUE)
  requireNamespace("iterators",quietly=TRUE)
  
  
  # Defines basic formulation for the model
  if(is.null(formula)){
    # formula=y~treatment+factor(time)+(1|cluster) ##categorical time effect
    formula=y~treatment+time+(1|cluster)##doesn't matter to use a linear term if timeeff=0
    if(design=="cohort") {formula <- update.formula(formula,.~.+(1|id))}
  }
  
  # Defines the name of the treatment variable (for which the main analysis is performed)
  if(exists("treat.name",where=exArgs)) {treatment <- exArgs$treat.name} else {treatment <- "treatment"}
  
  res <- list()
  tic <- proc.time()
  
  # For user-specified data generating processes
  ##deleted for convenience
  
  # For standard SWT data generating processes (uses make.swt)
  if(!exists("data",where=exArgs)) {
    if(exists("X",where=exArgs)) {
      X=exArgs$X
      row.names(X) <- sample(1:I,I)
      colnames(X) <- c("Baseline",paste0("Time ",1:J))
    } else {
      X=NULL
    }
    
    if (method=="lme") {
      signif <- vector()
      
      for (i in 1:n.sims) {
        
        fake.data <- make.swt(I=I,J=J,H=H,K=K,design=design,mu=mu,b.trt=b.trt,
                              b.time=b.time,sigma.y=sigma.y,sigma.e=sigma.e,
                              rho=rho,sigma.a=sigma.a,rho.ind=rho.ind,sigma.v=sigma.v,
                              X=X,family=family,natural.scale=natural.scale)
        
        # If the formula contains random effects, run lmer
        check.random.effect <- !is.null(findbars(formula))
        if(check.random.effect==TRUE) {
          if(family=="gaussian") {
            m <- lme4::lmer(formula, data=fake.data,
                            control = lmerControl(calc.derivs = FALSE, optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
          } else {
            m <- lme4::glmer(formula, data=fake.data,family=family,
                             control = glmerControl(calc.derivs = FALSE, optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
            
          }
          which.treat <- which(names(fixef(m))==treatment)
          theta <- lme4::fixef(m)[which.treat]
          sd.theta <- summary(m)$coefficients[which.treat,2]
          Vcov <- vcov(m, useScale = FALSE)
          se <- sqrt(diag(Vcov))
          betas <- lme4::fixef(m)
          tval <- betas / se
          ### NB: issue with p-values in random effect models (cannot determine df easily)
          ### An alternative (as recommended by Bates et al) is to use CIs instead of p-values
          ci.fixed <- confint(m,names(fixef(m)),method="Wald",level=1-sig.level)
          p0 <- mu
          p1 <- mu+ni.margin
         
            ##if upper limit of treatment effect is > log(OR), then inferior (not significant)
            ##if upper limit of treatment effect is <= log(OR), then non-inferior (significant)
            ##therefore, signif is as.numeric(<=)
            signif[i] <- as.numeric(ci.fixed[2,2]<log(p1*(1-p0)/((1-p1)*p0)))

          ### This estimates the p-value based on Normal approximation --- may not be too correct...
          ### See: http://mindingthebrain.blogspot.co.uk/2014/02/three-ways-to-get-parameter-specific-p.html
          pval <- 2*pnorm(abs(tval), lower.tail = FALSE)
          pvalue <- pval[which.treat] 
          ###signif <- pvalue < alpha
          VC <- as.data.frame(lme4::VarCorr(m))
          rnd.eff.sd <- VC[,which(colnames(VC)=="sdcor")]
          # Finds the rows that needs to be reported
          to.report <- c(which(is.na(VC[,3]==T), which(VC[,1]=="Residual")))
          rnd.eff.sd <- VC[to.report,which(colnames(VC)=="sdcor")]
          names(rnd.eff.sd)[which(is.na(VC[to.report,3])==T)] <- paste(VC[to.report,"grp"],VC[to.report,"var1"])
          names(rnd.eff.sd) <- gsub(" NA","",names(rnd.eff.sd))
          if(family=="gaussian") {method <- "lmer"} else {method <- "glmer"}
        }
      }
    }
  }
  
  toc <- proc.time(); time2run <- (toc-tic)[3]; names(time2run) <- "Time to run (secs)"
  
  ######
  # see: https://stackoverflow.com/questions/41372927/isincompletecon-error-when-knitting-pdf
  #closeAllConnections()
  ######
  power <- mean(signif)
  
  power
}