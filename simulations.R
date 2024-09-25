###functions used for simulations------------------------------------------
library(rms)

generate.data.s1 <- function(seed=1234, n=30, cluster.limit=c(10,10), mu=1, vu=1,vr=1, hlf=F){
  set.seed(seed)
  #cluster size
  if(hlf){
    n1 <- ceiling(n/2)
    size.cluster <- c(rep(cluster.limit[1], n1), rep(cluster.limit[2], n - n1))
  }
  else{
    if(cluster.limit[1] == cluster.limit[2]){
      size.cluster <- rep(cluster.limit[1], n)
    }else size.cluster <- replicate(n, sample(cluster.limit[1]:cluster.limit[2], 1))
  }
  #generate data for each cluster
  x <- list()
  cluster <- list()
  id <- list()
  #generate cluster means
  u <- rnorm(n, mean = mu, sd = sqrt(vu))
  for(i in 1:n){
    n.cluster <- size.cluster[i]
    r <- rnorm(n.cluster, mean = 0, sd = sqrt(vr))
    x[[i]] <- r+u[i] 
    cluster[[i]] <- rep(i, n.cluster)
    id[[i]] <- 1:n.cluster
  }
  dat <- data.frame("x"= unlist(x),
                    "cluster"=as.factor(unlist(cluster)),
                    "id"=as.factor(unlist(id)))
  dat$x <- exp(dat$x)
  return(dat)
}

phi2 <- function(c1,c2,r){
  mvtnorm::pmvnorm(upper = c(c1,c2), mean = c(0,0), sigma = matrix(c(1,r,r,1), ncol=2))
}
var1 <- function(theta, m, n, rho){
  c1 <- c2 <- qnorm(theta)
  (theta * (1-theta) + 2 * (n-1) * (phi2(c1,c2,(1+rho)/2) - theta^2) + (n-1)^2 * (phi2(c1,c2,rho) - theta^2) + (m * 2 - 2) * n * (phi2(c1,c2,1/2) - theta^2 + (n-1) * (phi2(c1,c2,rho/2) - theta^2))) / (m^2 * n^2)
}
varsrs <- function(theta, m, rho){
  c1 <- c2 <- qnorm(theta)
  (theta * (1-theta) + (2 * m - 2) * (phi2(c1,c2,1/2) - theta^2))/(m^2)
}
var0 <- function(m, n, rho){
  (1/4 + 2*(n-1)*asin((1+rho)/2)/(2*pi) + (n-1)^2*asin(rho)/(2*pi)+(2*m-2)*n*(1/12+(n-1)*asin(rho/2)/(2*pi)))/(m^2*n^2)
}


simPower_probit <- function(iter, m, g, r, beta, idx.r, idx.beta, lbs){
  sim_num <- 10
  v0 <- var1(1/2, m, g, r)
  theta <- exp(beta) * (exp(beta)-beta-1) / (exp(beta)-1)^2
  muy <- qnorm(theta) * sqrt(2)
  pval1 <- pval2 <- b_hat <- rep(NA, sim_num)
  ans <- cbind(pval1, pval2, b_hat)
  vr <- 1-r
  vu <- r
  trial <- (iter - 1) * sim_num
  filename <- paste('output/probit', "-",lbs,'cls', g, 'r', idx.r, 'b', idx.beta, '-',iter,'.RData', sep="")
  for(i in 1:sim_num){
    print(i)
    datx <- generate.data.s1(seed=i+trial, n = m, cluster.limit = c(g,g), mu=0, vu=vu, vr=vr)
    daty <- generate.data.s1(seed=i+trial+1000, n = m, cluster.limit = c(g,g), mu=muy, vu=vu, vr=vr)
    X <- datx$x
    Y <- daty$x
    cls <- c(paste(datx$cluster, "x", sep=""), paste(daty$cluster, "y", sep=""))
    #Wilcoxon
    R <- rank(c(X, Y))
    N <- 2 * m
    w_hat <- sum(R[1:length(X)]) - m * g * (g * N + 1) / 2
    Ri <- tapply(R, cls, sum)
    w_hat_var <- m ^ 2 / (N * (N - 1)) * sum((Ri - g * (1 + g * N) / 2)^2)
    ts <- w_hat/sqrt(w_hat_var)
    pval1[i] <- pnorm(abs(ts), lower.tail = F) * 2
    #PO
    dat <- data.frame(trt = rep(c("x","y"), each=length(X)),
                      val = c(X,Y), cluster = cls)
    mod <- orm(val ~ trt, data = dat, family='probit', x=T, y=T, maxit=35)
    mod_robust <- robcov(fit=mod, cluster = dat$cluster)
    pval2[i] <- anova(mod_robust)[2, 'P']
    b_hat[i] <- mod_robust$coefficients[length(mod$coefficients)]
    ans <- cbind(pval1, pval2, b_hat)
    save(ans, file=filename)
  }
  
}


simPower_unbalance <- function(iter, m, glimit, r, beta, idx.r, idx.beta, lbs, hlf=T){
  sim_num <- 10
  theta <- exp(beta) * (exp(beta)-beta-1) / (exp(beta)-1)^2
  muy <- qnorm(theta) * sqrt(2)
  pval1 <- pval3 <- b_hat <- rep(NA, sim_num)
  ans <- cbind(pval1, pval3, b_hat)
  vr <- 1-r
  vu <- r
  trial <- (iter - 1) * sim_num
  filename <- paste('output-ub/',lbs,'cls', glimit[1], '-', glimit[2], 'r', idx.r, 'b', idx.beta, '-',iter,'.RData', sep="")
  for(i in 1:sim_num){
    print(i)
    datx <- generate.data.s1(seed=i+trial, n = m, cluster.limit = glimit, mu=0, vu=vu, vr=vr, hlf=hlf)
    daty <- generate.data.s1(seed=i+trial+1000, n = m, cluster.limit = glimit, mu=muy, vu=vu, vr=vr, hlf=hlf)
    X <- datx$x
    Y <- daty$x
    cls <- c(paste(datx$cluster, "x", sep=""), paste(daty$cluster, "y", sep=""))
    #Wilcoxon
    R <- rank(c(X, Y))
    w_hat <- sum(R[1:length(X)]) - sum(R) / 2 #- m * g * (g * N + 1) / 2
    #m_g/N_g = 1/2, since m_g = n_g (Rosner et al., 2003)
    m1 <- n1 <- ceiling(m/2); m2 <- n2 <- m - m1
    N1 <- m1 + n1; N2 <- m2 + n2
    cls.size <- tapply(R, cls, length)
    Rg <- tapply(R, cls, I)
    #cls size = glimit[1] in both arms
    idx <- names(cls.size)[cls.size == glimit[1]]
    Rsg <- sum(unlist(Rg[idx]))
    v1 <- sum(unlist(lapply(Rg[idx], function(i) (sum(i) - Rsg / N1)^2)))
    v1 <- m1 * n1 / N1 / (N1 - 1) * v1 
    #cls size = glimit[2] in both arms
    idx <- names(cls.size)[cls.size == glimit[2]]
    Rsg <- sum(unlist(Rg[idx]))
    v2 <- sum(unlist(lapply(Rg[idx], function(i) (sum(i) - Rsg / N2)^2)))
    v2 <- m2 * n2 / N2 / (N2 - 1) * v2
    ts <- w_hat/sqrt(v1 + v2)
    pval1[i] <- pnorm(abs(ts), lower.tail = F) * 2
    #PO
    dat <- data.frame(trt = rep(c("x","y"), each=length(X)),
                      val = c(X,Y), cluster = cls)
    mod <- orm(val ~ trt, data = dat, family='logistic', x=T, y=T, maxit=35)
    mod_robust <- robcov(fit=mod, cluster = dat$cluster)
    pval3[i] <- anova(mod_robust)[2, 'P']
    b_hat[i] <- mod_robust$coefficients[length(mod$coefficients)]
    ans <- cbind(pval1, pval3, b_hat)
    save(ans, file=filename)
  }
}


rGumbelMin <- function(n, mu=0, sigma=1){
  #loglog
  u <- runif(n, min=0, max=1)
  x <- mu + sigma*log(-log(1-u))
  return(x)
}

generate.data.extreme <- function(seed=1234, n=30, cluster.limit=c(10,10), mu=1, vu=1, vr=1, t=c('loglog','cloglog')){
  set.seed(seed)
  #cluster size
  if(cluster.limit[1] == cluster.limit[2]){
    size.cluster <- rep(cluster.limit[1], n)
  }
  else size.cluster <- replicate(n, sample(cluster.limit[1]:cluster.limit[2], 1))
  #generate data for each cluster
  x <- list()
  cluster <- list()
  id <- list()
  #generate cluster means
  u <- rnorm(n, mean = mu, sd = sqrt(vu))
  for(i in 1:n){
    n.cluster <- size.cluster[i]
    r <- runif(n.cluster, min=0, max=1) 
    if(t == 'loglog'){
      x[[i]] <- u[i] + vr*log(-log(1-r))
    }
    else{
      x[[i]] <- u[i] + vr*(-log(-log(1-r)))
    }
    cluster[[i]] <- rep(i, n.cluster)
    id[[i]] <- 1:n.cluster
  }
  dat <- data.frame("x"= unlist(x),
                    "cluster"=as.factor(unlist(cluster)),
                    "id"=as.factor(unlist(id)))
  return(dat)
}

simPower_loglog <- function(iter, m, g, r, beta, idx.r, idx.beta, sigma_loglog, lbs){
  sim_num <- 10
  theta <- exp(beta) * (exp(beta)-beta-1) / (exp(beta)-1)^2
  pval3 <- b_hat <- rep(NA, sim_num)
  ans <- cbind(pval3, b_hat)
  vu <- 1.6449 * r / (1 - r)
  vr <- sigma_loglog
  trial <- (iter - 1) * sim_num
  filename <- paste('output/loglog-',lbs,'cls', g, 'r', idx.r, 'b', idx.beta, '-',iter,'.RData', sep="")
  for(i in 1:sim_num){
    print(i)
    muy <- qnorm(theta)*sqrt(1.6449 * 2 + vu * 2)
    datx <- generate.data.extreme(seed=i+trial, n = m, cluster.limit = c(g,g), mu=0, vu=vu, vr=vr,t='loglog')
    daty <- generate.data.extreme(seed=i+trial+1000, n = m, cluster.limit = c(g,g), mu=muy, vu=vu, vr=vr,t='loglog')
    X <- datx$x
    Y <- daty$x
    cls <- c(paste(datx$cluster, "x", sep=""), paste(daty$cluster, "y", sep=""))
    dat <- data.frame(trt = rep(c("x","y"), each=length(X)),
                      val = c(X,Y), cluster = cls)
    mod <- orm(val ~ trt, data = dat, family='loglog', x=T, y=T, maxit=35)
    mod_robust <- robcov(fit=mod, cluster = dat$cluster)
    pval3[i] <- anova(mod_robust)[2, 'P']
    b_hat[i] <- mod_robust$coefficients[length(mod$coefficients)]
    ans <- cbind(pval3, b_hat)
    save(ans, file=filename)
  }
}

rGumbelMax<- function(n, mu=0, sigma=1){
  u <- runif(n, min=0, max=1)
  x <- mu + sigma*(-log(-log(u)))
  return(x)
}

simPower_cloglog <- function(iter, m, g, r, beta, idx.r, idx.beta, sigma_loglog, lbs){
  sim_num <- 10
  theta <- exp(beta) * (exp(beta)-beta-1) / (exp(beta)-1)^2
  pval3 <- b_hat <- rep(NA, sim_num)
  ans <- cbind(pval3, b_hat)
  vu <- 1.6449 * r / (1 - r)
  vr <- sigma_loglog
  trial <- (iter - 1) * sim_num
  filename <- paste('output/POcloglog-',lbs,'cls', g, 'r', idx.r, 'b', idx.beta, '-',iter,'.RData', sep="")
  if(file.exists(filename)){
    load(filename)
    pval3 <- ans[,1]; b_hat <- ans[,2]
  }
  if(sum(is.na(ans[,1])) > 0){
    for(i in (sum(!is.na(ans[,1]))+1):sim_num){
      print(i)
      muy <- qnorm(theta)*sqrt(1.6449 * 2 + vu * 2)
      datx <- generate.data.extreme(seed=i+trial, n = m, cluster.limit = c(g,g), mu=0, vu=vu, vr=vr,t='cloglog')
      daty <- generate.data.extreme(seed=i+trial+1000, n = m, cluster.limit = c(g,g), mu=muy, vu=vu, vr=vr,t='cloglog')
      X <- datx$x
      Y <- daty$x
      cls <- c(paste(datx$cluster, "x", sep=""), paste(daty$cluster, "y", sep=""))
      dat <- data.frame(trt = rep(c("x","y"), each=length(X)),
                        val = c(X,Y), cluster = cls)
      mod <- orm(val ~ trt, data = dat, family='cloglog', x=T, y=T, maxit=35)
      mod_robust <- robcov(fit=mod, cluster = dat$cluster)
      pval3[i] <- anova(mod_robust)[2, 'P']
      b_hat[i] <- mod_robust$coefficients[length(mod$coefficients)]
      ans <- cbind(pval3, b_hat)
      save(ans, file=filename)
    }
  }
}

simFP <- function(iter, m, g, r, idx.r){
  sim_num <- 10
  v0 <- var1(1/2, m, g, r)
  pval1 <- pval2 <- pval3 <- b_hat <- rep(NA, sim_num)
  vr <- 1-r
  vu <- r
  trial <- (iter - 1) * sim_num
  filename <- paste('output/FP-n',m,'cls', g, 'r', idx.r, '-',iter,'.RData', sep="")
  for(i in 1:sim_num){
    print(i)
    datx <- generate.data.s1(seed=i+trial, n = m, cluster.limit = c(g,g), mu=0, vu=vu, vr=vr)
    daty <- generate.data.s1(seed=i+trial+1000, n = m, cluster.limit = c(g,g), mu=0, vu=vu, vr=vr)
    X <- datx$x
    Y <- daty$x
    cls <- c(paste(datx$cluster, "x", sep=""), paste(daty$cluster, "y", sep=""))
    #Wilcoxon
    R <- rank(c(X, Y))
    N <- 2 * m 
    w_hat <- sum(R[1:length(X)]) - m * g * (g * N + 1) / 2
    Ri <- tapply(R, cls, sum)
    w_hat_var <- m ^ 2 / (N * (N - 1)) * sum((Ri - g * (1 + g * N) / 2)^2)
    ts <- w_hat/sqrt(w_hat_var)
    pval1[i] <- pnorm(abs(ts), lower.tail = F) * 2
    #Wilcoxon theta
    theta_hat <- sum(apply(outer(X, Y, '-'), 1:2, fu)) / (m^2*g^2)
    tst <- (theta_hat - 1/2) / sqrt(v0)
    pval2[i] <- pnorm(abs(tst), lower.tail = F) * 2
    #PO
    dat <- data.frame(trt = rep(c("x","y"), each=length(X)),
                      val = c(X,Y), cluster = cls)
    mod <- orm(val ~ trt, data = dat, family='logistic', x=T, y=T, maxit=35)
    mod_robust <- robcov(fit=mod, cluster = dat$cluster)
    pval3[i] <- anova(mod_robust)[2, 'P']
    b_hat[i] <- mod_robust$coefficients[length(mod$coefficients)]
    ans <- cbind(pval1, pval3, pval2, b_hat)
    save(ans, file=filename)
  }
}

#empirically obtain true values of treatment effect and rank ICC of loglog/cloglog distributed data 
generate.data.extreme.emp <- function(seed=1234, n=30, cluster.limit=c(10,10), mu=1, vu=1, vr=1){
  set.seed(seed)
  #cluster size
  if(cluster.limit[1] == cluster.limit[2]){
    size.cluster <- rep(cluster.limit[1], n)
  }
  else size.cluster <- replicate(n, sample(cluster.limit[1]:cluster.limit[2], 1))
  #generate data for each cluster
  x.loglog <- x.cloglog <- list()
  cluster <- list()
  id <- list()
  #generate cluster means
  u <- rnorm(n, mean = mu, sd = sqrt(vu))
  for(i in 1:n){
    n.cluster <- size.cluster[i]
    r <- runif(n.cluster, min=0, max=1)  
    x.loglog[[i]] <- u[i] + vr*log(-log(1-r))
    x.cloglog[[i]] <- u[i] + vr*(-log(-log(1-r)))
    cluster[[i]] <- rep(i, n.cluster)
    id[[i]] <- 1:n.cluster
  }
  dat <- data.frame("xloglog"= unlist(x.loglog),
                    "xcloglog"= unlist(x.cloglog),
                    "cluster"=as.factor(unlist(cluster)),
                    "id"=as.factor(unlist(id)))
  return(dat)
}
betas <- c(0.5, 1, 1.5)
ans <- list()
for(idx.beta in seq_along(betas)){
  beta <- betas[idx.beta]
  theta <- exp(beta) * (exp(beta)-beta-1) / (exp(beta)-1)^2
  mu.x <- 0
  ri <- seq(0, 0.9, 0.1)
  us <- 1.6449 * ri / (1 - ri)
  sigma_extreme <- 1
  rs <- rep(sigma_extreme, length(ri))
  k <- length(us)
  g <- list()
  N <- 1000000
  for(i in seq_along(ri)){
    print(i)
    mu.y <- qnorm(theta)*sqrt(1.6449 * 2 + us[i] * 2)
    datx <- generate.data.extreme.emp(seed = i, n = N, cluster.limit = c(100, 100), mu = mu.x, vu = us[i], vr = rs[i])
    daty <- generate.data.extreme.emp(seed = i + length(ri), n = N, cluster.limit = c(100, 100), mu = mu.y, vu = us[i], vr = rs[i])
    ####repeated sampling clusters with replacement
    cls <- sample(unique(datx$cluster), size = N, replace = T)
    pairs <- list()
    lx <- tapply(1:nrow(datx), datx$cluster, I)
    for(j in seq_along(cls)){
      ix <- lx[[cls[j]]]
      pairs[[j]] <- sample(ix, 2, replace = F)
    }
    pairs <- do.call(rbind, pairs)
    
    idx1 <- sample(1:nrow(datx), size = N, replace = T)
    idx2 <- sample(1:nrow(datx), size = N, replace = T)
    
    est <- matrix(0, ncol = 3, nrow = 2)
    print("EST")
    #####################loglog
    xi <- datx$xloglog
    yi <- daty$xloglog
    ####rank ICC
    ex <- cor(xi[pairs[,1]], xi[pairs[,2]], method = "spearman")
    ey <- cor(yi[pairs[,1]], yi[pairs[,2]], method = "spearman")
    est[1, c(1:2)] <- c(ex, ey)
    #######theta
    #########empirical theta
    ki <- xi[idx1] - yi[idx2]
    est[1, 3] <- (sum(ki<0) + sum(ki==0)/2) / N
    
    #####################cloglog
    xi <- datx$xcloglog
    yi <- daty$xcloglog
    ####rank ICC
    ex <- cor(xi[pairs[,1]], xi[pairs[,2]], method = "spearman")
    ey <- cor(yi[pairs[,1]], yi[pairs[,2]], method = "spearman")
    est[2, c(1:2)] <- c(ex, ey)
    #######theta
    #########empirical theta
    ki <- xi[idx1] - yi[idx2]
    est[2, 3] <- (sum(ki<0) + sum(ki==0)/2) / N 
    
    g[[i]] <- est
    ans[[as.character(beta)]] <- g
    save(ans, file = "output/CRCT-extreme-values.RData")
  }
}

approxBeta <- function(t, tol=1e-6, maxIter=100){
  betas <- seq(-5,5,0.1)
  betas <- betas[betas!=0]
  theta <- exp(betas) * (exp(betas)-betas-1) / (exp(betas)-1)^2
  u0 <- min(betas[theta > t]) 
  l0 <- max(betas[theta < t])
  i <- 0; d <- 10
  while(i < maxIter & d > tol){
    b <- (u0+l0) / 2
    tb <- exp(b) * (exp(b)-b-1) / (exp(b)-1)^2
    if(tb > t) u0 <- b
    else l0 <- b
    d <- abs(tb - t)
    i <- i + 1
  }
  return(b)
}


simPower_ordinal <- function(iter, m, g, r, gr, beta, idx.r, idx.beta, lbs){
  ls <- c(3, 5, 10)
  sim_num <- 10
  theta <- exp(beta) * (exp(beta)-beta-1) / (exp(beta)-1)^2
  muy <- qnorm(theta) * sqrt(2)
  lans <- list()
  vr <- 1-r
  vu <- r
  trial <- (iter - 1) * sim_num
  filename = paste('output/qX-ordinal-', lbs, 'cls', g, 'r', idx.r, 'b', idx.beta, '-',iter,'.RData', sep="")
  for(i in 1:sim_num){
    pval2 <- b_hat <- rep(NA, length(ls))
    for(li in seq_along(ls)){
      mi <- m[li]
      datx <- generate.data.s1(seed=i+trial, n = mi, cluster.limit = c(g,g), mu=0, vu=vu, vr=vr)
      daty <- generate.data.s1(seed=i+trial+1000, n = mi, cluster.limit = c(g,g), mu=muy, vu=vu, vr=vr)
      
      cls <- c(paste(datx$cluster, "x", sep=""), paste(daty$cluster, "y", sep=""))
      l <- ls[li]
      xb <- qnorm(seq(1/l, 1-1/l, 1/l))
      bi <- c(-Inf, sort(xb), Inf)
      nls <- length(bi) - 1
      xi <- as.numeric(cut(datx$x, breaks = bi, labels = 1:nls))
      yi <- as.numeric(cut(daty$x, breaks = bi, labels = 1:nls))
      
      #wilcoxon
      R <- rank(c(xi, yi))
      N <- 2 * mi
      w_hat <- sum(R[1:length(xi)]) - mi * g * (g * N + 1) / 2
      Ri <- tapply(R, cls, sum)
      w_hat_var <- mi ^ 2 / (N * (N - 1)) * sum((Ri - g * (1 + g * N) / 2)^2)
      ts <- w_hat/sqrt(w_hat_var)
      pval1[li] <- pnorm(abs(ts), lower.tail = F) * 2
      #PO
      dat <- data.frame(trt = rep(c("x","y"), each=length(xi)),
                        val = c(xi,yi), cluster = cls)
      mod <- orm(val ~ trt, data = dat, family='logistic', x=T, y=T, maxit=35)
      mod_robust <- robcov(fit=mod, cluster = dat$cluster)
      pval3[li] <- anova(mod_robust)[2, 'P']
      b_hat[li] <- mod_robust$coefficients[length(mod$coefficients)]
      
      ans <- cbind(pval1, pval2, b_hat)
    }
    lans[[i]] <- ans
    save(lans, file=filename)
  }
}

#empirically obtain true values of theta and rank ICC of ordinal data 
betas <- c(0.5, 1, 1.5)
ans <- list()
for(idx.beta in seq_along(betas)){
  beta <- betas[idx.beta]
  theta <- exp(beta) * (exp(beta)-beta-1) / (exp(beta)-1)^2
  mu <- qnorm(theta)*sqrt(2)
  mu.x <- 0
  mu.y <- mu
  ri <- seq(0, 1, 0.1)
  us <- ri
  rs <- 1 - ri
  k <- length(us)
  g <- pt <- list()
  ls <- c(3, 5, 10)
  N <- 1000000
  for(i in seq_along(ri)){
    print(i)
    datx <- generate.data.s1(seed = i, n = N, cluster.limit = c(100, 100), mu = mu.x, vu = us[i], vr = rs[i])
    daty <- generate.data.s1(seed = i + length(ri), n = N, cluster.limit = c(100, 100), mu = mu.y, vu = us[i], vr = rs[i])
    ####repeated sampling clusters with replacement
    cls <- sample(unique(datx$cluster), size = N, replace = T)
    pairs <- list()
    lx <- tapply(1:nrow(datx), datx$cluster, I)
    for(j in seq_along(cls)){
      ix <- lx[[cls[j]]]
      pairs[[j]] <- sample(ix, 2, replace = F)
    }
    pairs <- do.call(rbind, pairs)
    
    idx1 <- sample(1:nrow(datx), size = N, replace = T)
    idx2 <- sample(1:nrow(datx), size = N, replace = T)
    # k2 <- datx$x[idx1] - daty$x[idx2]
    # pt[[i]] <- (sum(k2<0) + sum(k2==0)/2) / N
    
    # est <- matrix(0, ncol = 5, nrow = 3)
    est <- matrix(0, ncol = 3, nrow = 3)
    print("EST")
    for(li in seq_along(ls)){
      l <- ls[li]
      print(l)
      # xyb <- quantile(c(rnorm(1e6,  mu.x, 1), rnorm(1e6,  mu.y, 1)), seq(1/l, 1-1/l, 1/l))
      # xyb <- quantile(c(datx$x, daty$x), seq(1/l, 1-1/l, 1/l))
      # c(quantile(datx$x, seq(1/l, 1-1/l, 1/l)), quantile(daty$x, seq(1/l, 1-1/l, 1/l)))      
      # bi <- c(-Inf, sort(xyb), Inf)
      xb <- qnorm(seq(1/l, 1-1/l, 1/l))
      bi <- c(-Inf, sort(xb), Inf)
      nls <- length(bi) - 1
      xi <- as.numeric(cut(datx$x, breaks = bi, labels = 1:nls))  
      ex <- cor(xi[pairs[,1]], xi[pairs[,2]], method = "spearman")
      # ex2 <- cor(xi[pairs[,1]], xi[pairs[,2]], method = "pearson")
      #######y
      yi <- as.numeric(cut(daty$x, breaks = bi, labels = 1:nls))  
      ey <- cor(yi[pairs[,1]], yi[pairs[,2]], method = "spearman")
      # eys <- cor(yi[pairs[,1]], yi[pairs[,2]], method = "pearson")
      # est[li, c(1:2, 4:5)] <- c(ex, ey, ex2, ey2)
      est[li, c(1:2)] <- c(ex, ey)
      #######theta
      #########empirical theta
      ki <- xi[idx1] - yi[idx2]
      est[li, 3] <- (sum(ki<0) + sum(ki==0)/2) / N
    }
    g[[i]] <- est
    # ans[[as.character(beta)]] <- list(g, unlist(pt))
    ans[[as.character(beta)]] <- g
    save(ans, file = "output/CRCT-discrete.RData")
  }
}


###simulations for continuous data------------------------------------------
######probit link
###calculate sample sizes
a <- 0.05 #significance level
b <- 1 - 0.9 
lmg <- list()
betas <- c(0.1, 0.5, 1, 1.5, 2.5, 4) #log OR
rhos <- seq(0, 0.9, 0.1) #ICC
n <- 5 #number of clusters 
for(j in seq_along(betas)){
  beta <- betas[j]
  mg <- rep(NA, length(rhos))
  theta <- exp(beta) * (exp(beta)-beta-1) / (exp(beta)-1)^2
  for(i in seq_along(rhos)){
    rho <- rhos[i]
    ga <- 6 * asin(rho/2)/pi
    Dg <- 6 * (qnorm(1-a/2) + qnorm(1-b))^2 / beta^2 * (1 + (n - 1) * ga)
    mg[i] <- sqrt(1/4 + Dg^2/4) + Dg/2
  }
  lmg[[j]] <- mg
}

###Power
for(iter in 1:100){
  n <- 5 #number of clusters
  for(j in seq_along(betas)){
    beta <- betas[j]
    for(i in seq_along(rhos)){
      print(paste("j:",j," ,i:",i, sep=""))
      m <- ceiling(lmg[[j]][i] / n)
      rho <- rhos[i]
      simPower_probit(iter, m, n, rho, beta, i, j, "") 
    }
  }
}

######loglog and cloglog link
###calculate sample sizes 
load("CRCT-extreme-values.RData")
betas <- c(0.5, 1, 1.5)
rhos <- seq(0, 0.9, 0.1)
thetas <- betas.logit <- list()
gIs <- list()
for(i in seq_along(betas)){
  ansi <- ans[[i]]
  gIs[[i]] <- lapply(ansi, function(x){
    t <- rowMeans(x[,1:2])
    t[t < 0] <- 0
    t
  })
  thetas[[i]] <- do.call(rbind, lapply(ansi, function(x) x[,3]))
  betas.logit[[i]] <- apply(thetas[[i]], c(1,2), function(x) approxBeta(x))
}

a <- 0.05
b <- 1 - 0.9
lm_loglog <- lm_cloglog <- list()
betas <- c(0.5, 1, 1.5)
rhos <- seq(0, 0.9, 0.1)
n <- 5
for(j in seq_along(betas)){
  betai.logit <- betas.logit[[j]]
  m.loglog <- m.cloglog <- rep(NA, length(rhos))
  gI <- gIs[[j]]
  for(i in seq_along(rhos)){
    beta.logit <- betai.logit[i,]
    gIi <- gI[[i]]
    Dg.loglog <- 6 * (qnorm(1-a/2) + qnorm(1-b))^2 / beta.logit[1]^2 * (1 + (n - 1) * gIi[1])
    m.loglog[i] <- sqrt(1/4 + Dg.loglog^2/4) + Dg.loglog/2
    Dg.cloglog <- 6 * (qnorm(1-a/2) + qnorm(1-b))^2 / beta.logit[2]^2 * (1 + (n - 1) * gIi[2])
    m.cloglog[i] <- sqrt(1/4 + Dg.cloglog^2/4) + Dg.cloglog/2
  }
  lm_loglog[[j]] <- ceiling(m.loglog / n)
  lm_cloglog[[j]] <- ceiling(m.cloglog / n)
}

ms <- list(lm_loglog, lm_cloglog, gIs, thetas, betas.logit)
save(ms, file = "extreme_values_size_cls5.Rdata")

###power 
load("extreme_values_size_cls5.Rdata")
ms.loglog <- ms[[1]]
ms.cloglog <- ms[[2]]
for(iter in 1:100){
  betas <- c(0.1, 0.5, 1, 1.5, 2.5, 4)
  rhos <- seq(0, 0.9, 0.1)
  n <- 5
  sigma_loglog <- 1
  for(j in seq_along(betas)[3:4]){
    beta <- betas[j]
    for(i in seq_along(rhos)){
      rho <- rhos[i]
      print(paste("j:",j," ,i:",i, sep=""))
      m.loglog <- ms.loglog[[j-1]][i]
      print(paste("loglog", m.loglog, sep=": "))
      simPower_loglog(iter, m.loglog, n, rho, beta, i, j, sigma_loglog, "")
      m.cloglog <- ms.cloglog[[j-1]][i]
      print(paste("cloglog", m.cloglog, sep=": "))
      simPower_cloglog(iter, m.loglog, n, rho, beta, i, j, sigma_loglog, "")
    }
  }
}


#######unbalance data
#sample size
a <- 0.05 #significance level
b <- 1 - 0.9 
lmg <- list()
betas <- c(0.1, 0.5, 1, 1.5, 2.5, 4) #log OR
rhos <- seq(0, 0.9, 0.1) #ICC
n <- 5 #number of clusters 
for(j in seq_along(betas)){
  beta <- betas[j]
  mg <- rep(NA, length(rhos))
  theta <- exp(beta) * (exp(beta)-beta-1) / (exp(beta)-1)^2
  for(i in seq_along(rhos)){
    rho <- rhos[i]
    ga <- 6 * asin(rho/2)/pi
    Dg <- 6 * (qnorm(1-a/2) + qnorm(1-b))^2 / beta^2 * (1 + (n - 1) * ga)
    mg[i] <- sqrt(1/4 + Dg^2/4) + Dg/2
  }
  lmg[[j]] <- mg
}
#power
for(iter in 1:100){
  for(j in seq_along(betas)){
    beta <- betas[j]
    for(i in seq_along(rhos)){
      print(paste("j:",j," ,i:",i, sep=""))
      m <- ceiling(lmg[[j]][i] / n)
      rho <- rhos[i]
      simPower_unbalance(iter, m, c(10, 30), rho, beta, i, j, "", T)
      simPower_unbalance(iter, m, c(15, 25), rho, beta, i, j, "", F)
    }
  }
}

#########R&G sample sizes 
a <- 0.05
b <- 1 - 0.9
lmrg_r <- list()
betas <- c(0.1, 0.5, 1, 1.5, 2.5, 4)#log OR
rhos <- seq(0, 0.9, 0.1) #icc
n <- 5 #number of clusters
for(j in seq_along(betas)){
  beta <- betas[j]
  mrg_r <- rep(NA, length(rhos))
  theta <- exp(beta) * (exp(beta)-beta-1) / (exp(beta)-1)^2
  for(i in seq_along(rhos)){
    rho <- rhos[i]
    mrg_r[i] <- approxSize(n, theta, a, 1 - b, rho) * n 
  }
  lmrg_r[[j]] <- mrg_r 
}

####R&G ower
for(iter in 1:100){
  n <- 5
  for(j in seq_along(betas)){
    beta <- betas[j]
    for(i in seq_along(rhos)){
      print(paste("j:",j," ,i:",i, sep=""))
      m <- ceiling(lmrg_r[[j]][i] / n)
      rho <- rhos[i]
      simPower(iter, m, n, rho, beta, i, j, "RG-")
    }
  }
  
}

####R&G under unbalance
for(iter in 1:100){
  for(j in seq_along(betas)){
    beta <- betas[j]
    for(i in seq_along(rhos)){
      print(paste("j:",j," ,i:",i, sep=""))
      m <- ceiling(lmrg_r[[j]][i] / n)
      rho <- rhos[i]
      simPower_unbalance(iter, m, c(15, 25), rho, beta, i, j, "RG-")
      simPower_unbalance(iter, m, c(5, 35), rho, beta, i, j, "RG-")
    }
  }
  
}


###simulations for ordinal data------------------------------------------
#sample sizes 
betas <- c(0.5, 1, 1.5)
ls <- c(3,5,10)
nl3 <- matrix(NA, ncol = length(ls), nrow = length(betas))
for(idx.beta in seq_along(betas)){
  beta <- betas[idx.beta]
  theta <- exp(beta) * (exp(beta)-beta-1) / (exp(beta)-1)^2
  mu <- qnorm(theta) * sqrt(2)
  for(idx.l in seq_along(ls)){
    l <- ls[idx.l]
    xb <- qnorm(seq(1/l, 1-1/l, 1/l))
    bi <- c(-Inf, sort(xb), Inf)
    xn <- sapply(1:(length(bi)-1), function(i) pnorm(bi[i+1])-pnorm(bi[i]))
    yn <- sapply(1:(length(bi)-1), function(i) pnorm(bi[i+1], mu)-pnorm(bi[i], mu))
    ln <- (xn + yn)/2
    nl3[idx.beta, idx.l] <- 1 - sum(ln^3)
  }
}

approxBeta <- function(t, tol=1e-6, maxIter=100){
  betas <- seq(-5,5,0.1)
  betas <- betas[betas!=0]
  theta <- exp(betas) * (exp(betas)-betas-1) / (exp(betas)-1)^2
  u0 <- min(betas[theta > t])
  l0 <- max(betas[theta < t])
  i <- 0; d <- 10
  while(i < maxIter & d > tol){
    b <- (u0+l0) / 2
    tb <- exp(b) * (exp(b)-b-1) / (exp(b)-1)^2
    if(tb > t) u0 <- b
    else l0 <- b
    d <- abs(tb - t)
    i <- i + 1
  }
  return(b)
}

ordinalSize <- function(a, b, lor, nls, r, g){
  6 * (qnorm(1-a/2) + qnorm(1-b))^2 * (1 + (g - 1) * r) / lor^2 / nls
}

load("accre/CRCT-discrete.RData")
r <- list()
t <- matrix(NA, ncol = 3, nrow = length(ans))
for(i in seq_along(ans)){
  ans.sub <- ans[[i]]
  r[[i]] <- do.call(rbind, lapply(ans.sub, function(x) rowMeans(x[,1:2])))
  t[i, ] <- colMeans(do.call(rbind, lapply(ans.sub, function(x) x[,3])))
}
betas.ord <- apply(t, c(1,2), function(x) approxBeta(x))


a <- 0.05
b <- 1 - 0.9
g <- 5
ms <- list()
for(i in seq_along(r)){
  beta <- betas.ord[i,]
  theta <- exp(beta) * (exp(beta)-beta-1) / (exp(beta)-1)^2
  rhos <- r[[i]]
  mi <- mrgi <- matrix(NA, ncol = 3, nrow = nrow(rhos))
  for(j in 1:nrow(rhos)){
    mi[j, 1] <- ceiling(ordinalSize(a, b, beta[1], nl3[i, 1], rhos[j,1], g) / g)
    mi[j, 2] <- ceiling(ordinalSize(a, b, beta[2], nl3[i, 2], rhos[j,2], g) / g)
    mi[j, 3] <- ceiling(ordinalSize(a, b, beta[3], nl3[i, 3], rhos[j,3], g) / g)
  }
  ms[[i]] <- mi
}
msrs <- list(ms, r)
save(msrs, file = "ordinal_size_cls5.Rdata")


load("ordinal_size_cls5.Rdata")
for(iter in 1:100){
  ms <- msrs[[1]]
  rs <- msrs[[2]]
  a <- 0.05
  b <- 1 - 0.9
  betas <- c(0.1, 0.5, 1, 1.5, 2.5, 4)
  rhos <- seq(0, 0.9, 0.1)
  n <- 5 #number of clusters
  for(j in seq_along(betas)[2:4]){
    beta <- betas[j]
    for(i in seq_along(rhos)){
      print(paste("j:",j," ,i:",i, sep=""))
      m <- ceiling(ms[[j-1]][i,])
      gr <- rs[[j-1]][i,]
      rho <- rhos[i]
      simPower_ordinal(iter, m, n, rho, gr, beta, i, j, "")
    }
  }
}

