### code taken from http://www.bios.unc.edu/research/genomic_software/Multi-Tissue-eQTL/code3.html
### Set estimation parameters
maxIterations = 100;

### Load big matrix of z-scores
load(file='z-score.matrix.Rdata');
dim(big.matrix);

### Initialize parameters
{
  param = list();
  
  ### K - the number of tissues
  K = ncol(big.matrix);
  
  ### Delta - null covariance matrix across tissues
  param$Delta = matrix(0.05, K, K);
  diag(param$Delta) = 1;

  ### Sigma - signal covariance matrix across tissues
  param$Sigma = matrix(3, K, K) + diag(K);
  
  ### P - the vector of probabilities
  param$P = rep(1/2^K, 2^K);
  
  ### Psubs - the vector of active eQTLs
  ###         for each element of P
  Psubs = vector('list',2^K);
  for( i in 1:2^K) {
    a<-2^((K-1):0);
    b<-2*a;
    Psubs[[i]] = as.double(((i-1) %% b)>=a);
  }
  rm(a,b,i);
  param$Psubs = Psubs;
  rm(Psubs);

  ### loglik - the initial likelihood
  param$loglik = -Inf;
  rm(K);
}

### The function does a single iteration of the estimation procedure
DoIteration = function(big.matrix, param) {
  ### extract current model parameters
  K = ncol(big.matrix);
  m = nrow(big.matrix);
  Delta = param$Delta;
  Sigma = param$Sigma;
  P = param$P;
  Psubs = param$Psubs;

  ### The function for matrix power
  mat.power = function(mat,pow) {
    e = eigen(mat);
    V = e$vectors;
    return( V %*% diag(e$values^pow) %*% t(V) );
  }
    
  ### Start the timer
  tic = proc.time();
  
  ### variables to accumulate
  ### loglik - likelihood
  ### newP - marginal probabilities
  ### newDelta - the new Delta matrix
  ### newSigmaPlusDelta - Delta+Sigma
  cum.loglik = 0;
  cum.newP = 0;
  cum.newDelta = 0;
  cum.newSigmaPlusDelta = 0;
  
  ### Do calculations in slices of 10000 gene-SNP pairs
  step1 = 1e5L;
  for( j in 1:ceiling(m/step1)) {
    fr = step1*(j-1)+1;
    to = min(step1*j,m);
    X = big.matrix[fr:to, , drop=FALSE];
    
    ### likelihood for the slice
    prob = matrix(0, nrow(X), length(P));
    for( i in 1:length(Psubs)) {
      sigma_star = Delta + Sigma * tcrossprod(Psubs[[i]]);
      sigma_hfiv = mat.power(sigma_star,-0.5);
      sigma_dethfiv = (det(sigma_star))^(-0.5);
      w = (1/(2*pi)^(K/2)) * (P[i]*sigma_dethfiv);
      prob[,i] = exp(log(w) - colSums(tcrossprod(sigma_hfiv/sqrt(2),X)^2) ) ;
    }
    
    cum.loglik = cum.loglik + sum(log(rowSums(prob)));
    
    ### Normalize probabilities for each gene-SNP pair
    ### to add up to 1
    prob = prob / rowSums(prob);
    ### new vector of P - tissue specificity probabilities
    cum.newP = cum.newP + colSums(prob);
    
    cum.newDelta = cum.newDelta +
                    crossprod(X*sqrt(prob[,1]));
    cum.newSigmaPlusDelta = cum.newSigmaPlusDelta +
                    crossprod(X*sqrt(prob[,length(P)]));
  }
  
  {
    ### Calculate Delta from the cumulative sum
    Delta = cum.newDelta / cum.newP[1];
    ### normalize to force the diagonal to 1
    Delta = Delta * tcrossprod(sqrt(1 / diag(Delta)));
    
    ### Same with Sigma
    Sigma = cum.newSigmaPlusDelta / tail(cum.newP,1) - Delta;
    e = eigen(Sigma);
    if( any(e$values<0) ) {
      Sigma = e$vectors %*% diag(pmax(e$values,0)) %*% t(e$vectors);
    }  
  }
  P = cum.newP / sum(cum.newP);
  
  toc = proc.time();
  return(list(Delta = Delta, Sigma = Sigma, P = P,
              Psubs = Psubs, loglik = cum.loglik, time = toc-tic));
}

### The 'paralist' list vector
### will store model estimates at each iteration
paralist = vector('list',maxIterations+1);
paralist[[1]] = param;
rm(param);

### Perform up to 'maxIterations' iteration
for( i in 2:length(paralist) ) {
  paralist[[i]] = DoIteration(big.matrix=big.matrix, param=paralist[[i-1]]);
  cat(i, '\t', paralist[[i]]$loglik-paralist[[i-1]]$loglik,
         '\t', paralist[[i]]$time[3],'\n');
  if(i>10)
    if(paralist[[i]]$loglik < paralist[[i-1]]$loglik)
      break;
}
paralist = paralist[!sapply(paralist, is.null)];

### Save the results
save(list='paralist', file='paralist.Rdata');

