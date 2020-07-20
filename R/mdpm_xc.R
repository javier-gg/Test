mdpm_xc<-function(y, L, prior, nsim, nburn, nthin = 1, start = NULL, scale = TRUE){
  ## Convert data into a matrix
  y<-cbind(y)

  ## Auxiliar variables
  n<-nrow(y)
  p<-ncol(y)

  ## Check prior information
  if("mu_0" %in% names(prior)){
    ## Extract hyperparameters
    mu_0<-prior$mu_0
  }
  else{
    if(all(c("m", "S") %in% names(prior))){
      ## Initial value
      mu_0<-rep(0, p)

      ## Hyperparameters
      m<-prior$m
      S<-prior$S
      S_1<-solve(S)
    }
    else{
      stop("If mu_0 is not supplied, hyperparameters m and S must be provided.")
    }
  }
  if("V_0" %in% names(prior)){
    ## Extract hyperparameters
    V_0<-prior$V_0
    V_0_1<-solve(V_0)
    v0<-FALSE
  }
  else{
    if(all(c("r", "Q") %in% names(prior))){
      ## Initial value
      V_0<-diag(p)
      v0<-TRUE

      ## Hyperparameters
      r<-prior$r
      Q<-prior$Q
    }
    else{
      stop("If V_0 is not supplied, hyperparameters r and Q must be provided.")
    }
  }
  if("nu_0" %in% names(prior)){
    ## Extract hyperparameters
    nu_0<-prior$nu_0
  }
  else{
    stop("Hyperparameter nu_0 must be provided.")
  }
  if("Psi_0" %in% names(prior)){
    ## Extract hyperparameters
    Psi_0<-prior$Psi_0
  }
  else{
    if(all(c("nu_1", "Psi_1") %in% names(prior))){
      ## Initial value
      Psi_0<-diag(p)

      ## Hyperparameters
      nu_1<-prior$nu_1
      Psi_1<-prior$Psi_1
      Psi_1_1<-solve(Psi_1)
    }
    else{
      stop("If Psi_0 is not supplied, hyperparameters nu_1 and Psi_1 must be provided.")
    }
  }
  if("alpha" %in% names(prior)){
    ## Extract hyperparameters
    alpha<-prior$alpha
  }
  else{
    if(all(c("a_alpha", "b_alpha") %in% names(prior))){
      ## Initial value
      alpha<-1

      ## Hyperparameters
      a_alpha<-prior$a_alpha
      b_alpha<-prior$b_alpha
    }
    else{
      stop("If alpha is not supplied, hyperparameters a_alpha and b_alpha must be provided.")
    }
  }

  ## Output objects
  params<-list()

  ## Scale the data
  y0<-y
  if(scale){
    m0<-colMeans(y0)
    ev<-eigen(cov(y))
    B<-ev$vectors%*%
      diag(ev$values, nrow = length(ev$values))^(0.5)%*%
      t(ev$vectors)
    B_1<-solve(B)
    y<-t(matrix(apply(y, 1, function(x) B_1%*%(x - m0)), ncol = n))
  }

  ## Initial values
  v<-rep(1/L, L)
  v[L]<-1
  comp<-comp_aux<-list()
  if(is.null(start)){
    for(l in 1:L){
      comp[[l]]<-list(omega = 1/L, mu = colMeans(y),
                      sigma = cov(y))
    }
  }
  else{
    if(all(c("omega", "mu", "sigma") %in% names(start))){
      if(sum(omega) == 1){
        omega<-start$omega
      }
      else{
        stop("Weights must sum up to one.")
      }
      mu<-as.matrix(start$mu, ncol = L)
      sigma<-start$sigma
      for(l in 1:L){
        comp[[l]]<-list(omega = omega[l], mu = mu[, l],
                        sigma = sigma[[l]])
      }
    }
    else{
      stop("Initial values for omega, mu and sigma must be specified.")
    }
  }

  ## Auxiliar variables
  pi_z<-matrix(, nrow = n, ncol = L)

  for(s in 1:nsim){
    ## Compute the weights
    comp[[1]]$omega<-v[1]
    for(l in 2:L){
      comp[[l]]$omega<-v[l]*(1 - v[l - 1])*
        comp[[l - 1]]$omega/v[l - 1]
    }

    ## Sample z
    for(l in 1:L){
      pi_z[, l]<-comp[[l]]$omega*mvnfast::dmvn(y, comp[[l]]$mu,
                                               comp[[l]]$sigma)
    }
    z<-apply(pi_z, 1, function(x) sample(1:L, size = 1, prob = x))

    ## Number of observations per cluster
    nl<-tabulate(z, L)

    ## Sample v
    for(l in 1:(L - 1)){
      v[l]<-rbeta(1, nl[l] + 1, alpha + sum(nl[(l + 1):L]))
    }

    ## Sample mu and sigma
    if(v0) V_0_1<-solve(V_0)
    for(l in 1:L){
      yz<-matrix(y[z == l, ], ncol = p)

      ## Sample mu
      syz<-colSums(yz)
      sigma_1<-solve(comp[[l]]$sigma)
      VAR<-solve(V_0_1 + nl[l]*sigma_1)
      MN<-V_0_1%*%mu_0 + sigma_1%*%syz
      comp[[l]]$mu<-MASS::mvrnorm(1, VAR%*%MN, VAR)

      ## Sample sigma
      y_star<-sweep(yz, 2, comp[[l]]$mu, "-")
      yty<-t(y_star)%*%y_star
      comp[[l]]$sigma<-CholWishart::rInvWishart(1, nu_0 + nl[l],
                                                Psi_0 + yty)[, , 1]
    }

    ## Sample mu_0
    if(!("mu_0" %in% names(prior))){
      VAR0<-solve(S_1 + L*V_0_1)
      MN0<-S_1%*%m + V_0_1%*%rowSums(matrix(sapply(comp,
                                                   function(x) x$mu), nrow = p))
      mu_0<-MASS::mvrnorm(1, VAR0%*%MN0, VAR0)
    }

    ## Sample V_0
    if(!("V_0" %in% names(prior))){
      mu_mu0<-matrix(sapply(comp,
                            function(x) x$mu - mu_0), nrow = p)
      V_0<-CholWishart::rInvWishart(1, r + L, Q + mu_mu0%*%t(mu_mu0))[, , 1]
    }

    ## Sample Psi_0
    if(!("Psi_0" %in% names(prior))){
      ss<-Reduce("+", lapply(comp, function(x) solve(x$sigma)))
      Psi_0<-rWishart(1, nu_1 + L*nu_0,
                                   solve(Psi_1_1 + ss))[, , 1]
    }

    ## Sample alpha
    if(!("alpha" %in% names(prior))){
      alpha<-rgamma(1, a_alpha + L - 1,
                    b_alpha - sum(log(1 - v[1:(L - 1)])))
    }

    ## Re-scale parameters
    if(scale){
      for(l in 1:L){
        comp_aux[[l]]<-list(omega = comp[[l]]$omega,
                            mu = B%*%comp[[l]]$mu + m0,
                            sigma = B%*%comp[[l]]$sigma%*%t(B))
      }
    }
    else{
      comp_aux<-comp
    }

    ## Store results
    params[[s]]<-comp_aux
  }

  ## Burn-in
  params<-lapply(seq(nburn + 1, nsim, nthin), function(x) params[[x]])

  return(list(y = y0, params = params))
}
