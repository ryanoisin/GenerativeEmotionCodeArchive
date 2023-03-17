# jonashaslbeck@gmail.com; Ovt 21

# Function to get phi and psi matrices out of mlVAR

f_getPars <- function(out, names) {
  
  out_sum <- summary(out)
  p <- length(names)
  
  # Lagged effects
  phi <- matrix(NA, p, p)
  n_edge <- p^2
  for(i in 1:n_edge) {
    ind1 <- which(names==out_sum$temporal[i, 1])
    ind2 <- which(names==out_sum$temporal[i, 2])
    phi[ind1, ind2] <- out_sum$temporal[i, 4]
  }
  
  # Residual Pcors
  psi <- matrix(0, p, p)
  n_edge <- p*(p-1)/2
  for(i in 1:n_edge) {
    ind1 <- which(names==out_sum$contemporaneous[i, 1])
    ind2 <- which(names==out_sum$contemporaneous[i, 2])
    psi[ind1, ind2] <- out_sum$contemporaneous[i, 5]
  }
  psi <- psi + t(psi)
  diag(psi) <- NA
  
  outlist <- list("phi" = phi, 
                  "psi" = psi)
  
  return(outlist)
  
} # eoF

# function to fit n=1 VAR model to data
fitVAR <- function(data, roundpoints = 4, summaries = FALSE, changescore = FALSE) {
  
  p <- ncol(data)
  Phi <- matrix(NA, p, p)
  alpha <- rep(NA, p)
  residuals <- matrix(NA, nrow(data)-1, p)
  summary_list<- list()
  for(i in 1:p) {
    
    # Estimates the change-score version of the var model instead
    if(changescore == TRUE){
      y = data[-1,i] - data[-nrow(data), i]
    } else {y <- data[-1, i]}
    
    
    X <- data[-nrow(data), ]
    
    coefs <- lm(y ~ X)
    Phi[i, ] <- coefs$coefficients[-1]
    alpha[i] <- coefs$coefficients[1]
    residuals[,i] <- coefs$residuals
    if(summaries) summary_list[[i]] <- summary(coefs)
  }
  Psi <- cov(residuals)
  mu <- as.vector(solve(diag(p)-Phi)%*%alpha)
  coefs <- list(round(alpha, roundpoints),round(mu,roundpoints), round(Phi, roundpoints), round(Psi,roundpoints))
  names(coefs) <- c("intercepts", "means","phi","psi")
  if(summaries){ return(summary_list) } else { return(coefs) }
  
}


# generate random transition matrix
gen_tmat <- function(p = 5, type = "any"){
  mout <- matrix(NA, p,p)
  if(type == "any"){
    for(j in 1:p){
      q <- runif(p)
      q <- q/sum(q)
      mout[j,] <- q
    }
  }else if(type == "diagheavy"){
    P <- matrix(runif(p*p), p, p)
    diag(P) <- (colSums(P) + rowSums(P)) /(p-1)
    round(P, 2)
    
    Pn <- apply(P, 1, function(x) x/sum(x))
    mout <- t(round(Pn, 2))
  }
  mout
}

# function check whether phenomena are present or not

check_phenom <- function(phi, om){
  p <- nrow(phi)
  # Phenomena 1: autocorrelations are positive
  p1 <- sum(diag(phi) > 0)/p
  
  # Phenomena 2:
  # AR bigger in absolute value than cross-lag
  # we use "majority rule"
  cl <- phi[row(phi)!=col(phi)]
  ar <- diag(phi)
  
  # output proportion of al > cl comparisons which are true
  p2 <- sum(sapply(ar, function(i) abs(i) > abs(cl) ))/(length(cl)*length(ar))
  
  # Phenomena 3.1: Within valence positive effects
  lwith <- as.matrix(bdiag(matrix(TRUE,3,3), matrix(TRUE,3,3)))
  diag(lwith) <- FALSE # don't count diagonals
  
  # within positive
  p3a <- sum(phi[lwith] > 0) / 12
  
  # Phenomena 3.2:  Between valence effects negative
  lbet <- !as.matrix(bdiag(matrix(TRUE,3,3), matrix(TRUE,3,3)))
  p3b <- sum(phi[lbet] < 0) / 18
  
  # Phenomana 4.1: Residuals within positive
  p4a <- sum(om[lwith] > 0) / 12
  
  # Phenomana 4.2: Residuals between negative
  p4b <-  sum(om[lbet] < 0) / 18
  
  c(p1, p2, p3a, p3b, p4a, p4b)
}
