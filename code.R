### preprocess
metadata <- read.csv("~/Metadata.csv")
C8p <- read.csv("~/C8p_NN_MN.csv")
C18n <- read.csv("~/C18n_NN_MN.csv")
HILp <- read.csv("~/HILp_NN_MN.csv")
df <- rbind(C18n,C8p,HILp)
df[is.na(df)] <- 0
del_row <- which(unlist(sapply(1:nrow(df), function(c) length(which(df[c,]==0))>0))==TRUE)
df <- df[-c(del_row),]
name_CD <- metadata$Tube.A..Metabolomics[which(metadata$diagnosis=="CD")]
name_UC <- metadata$Tube.A..Metabolomics[which(metadata$diagnosis=="UC")]
name_nonIBD <- metadata$Tube.A..Metabolomics[which(metadata$diagnosis=="nonIBD")]
df_CD <- df[,unlist(sapply(1:length(name_CD), function(c)which(name_CD[c]==colnames(df))))]
df_UC <- df[,unlist(sapply(1:length(name_UC), function(c)which(name_UC[c]==colnames(df))))]
df_nonIBD <- df[,unlist(sapply(1:length(name_nonIBD), function(c)which(name_nonIBD[c]==colnames(df))))]



### power fit
power_equation <- function(x, power_par){ t(sapply(1:nrow(power_par),function(c) power_par[c,1]*x^power_par[c,2] ) )}
power_equation_base <- function(x, y){
  x <- as.numeric(x)
  y <- as.numeric(y)
  min_value = min(y[y!=0])
  
  lmFit <- lm( log( y + runif(1, min = 0, max = min_value))  ~ log(x))
  coefs <- coef(lmFit)
  a <- exp(coefs[1])
  b <- coefs[2]
  
  model <- try(nls(y~a*x^b,start = list(a = a, b = b),
                   control = nls.control(maxiter = 1e3, minFactor = 1e-200)))
  if( 'try-error' %in% class(model)) {
    result = NULL
  }
  else{
    result = model
  }
  return(result)
}
power_equation_all <- function(x,y, maxit=1e2){
  result <- power_equation_base(x,y)
  iter <- 1
  while( is.null(result) && iter <= maxit) {
    iter <- iter + 1
    try(result <- power_equation_base(x,y))
  }
  return(result)
}
power_equation_fit <- function(data, trans = log10, thread = 1) {
  data = data[,order(colSums(data))]
  if ( is.null(trans)) {  
    X = colSums(data)
    trans_data = data
  } else{
    X = trans(colSums(data)+1)
    trans_data = trans(data+1)
  }
  colnames(trans_data) = X
  
  core.number <- thread
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterExport(cl, c("power_equation_all", "power_equation_base", "trans_data", "X"), envir = environment())
  all_model = parLapply(cl = cl, 1:nrow(data), function(c) power_equation_all(X, trans_data[c,]))
  stopCluster(cl)
  
  
  names(all_model) = rownames(data)
  no = which(sapply(all_model, length)>=1)
  all_model2 = all_model[no]
  data2 = data[no,]
  trans_data2 = trans_data[no,]
  
  new_x = X
  power_par = t(vapply(all_model2, coef, FUN.VALUE = numeric(2), USE.NAMES = TRUE))
  power_fit = t(vapply(all_model2, predict, newdata = data.frame(x=new_x),
                       FUN.VALUE = numeric(length(new_x)), USE.NAMES = TRUE))
  
  colnames(power_fit) = new_x
  result = list(original_data = data2, trans_data = trans_data2,
                power_par = power_par, power_fit = power_fit,
                Time = X)
  return(result)
}


### cluster
power_equation <- function(par,x){
  y <- par[1]*x^par[2]
  y
}
power_par <- function(y,times){
  
  x <- as.numeric(times)
  y <- as.numeric(y)
  lmFit <- lm( log( y + runif(1, min = 0, max = 0.1))  ~ log(x))
  coefs <- coef(lmFit)
  a <- exp(coefs[1])
  b <- coefs[2]
  
  tmp <- c(a,b)
  par_est <- function(par,x){
    sum( (y - power_equation(par,x))^2 )
  }
  r <- optim(tmp,par_est,x = as.numeric(times),method = "Nelder-Mead")
  return(r$par)
}
get_init_par <- function(data,k){

  #get initial pars based on k-means
  init_cluster <- kmeans(data,centers = k,iter.max = 1000)
  cuM <- init_cluster$centers 
  cusd <- diag(cov(df))
  init_curve_par <- cbind(t(sapply(1:k, function(c) power_par(cuM[c,1:265],times = meta_index1))),
                          t(sapply(1:k, function(c) power_par(cuM[c,266:411],times = meta_index2))),
                          t(sapply(1:k, function(c) power_par(cuM[c,412:546],times = meta_index3))))
  
  init_SAD_par <- c(mean(cusd[1:265]),0.4,mean(cusd[266:411]),0.4,mean(cusd[412:546]),0.4)
  
  init_pro <- table(init_cluster$cluster)/nrow(data)
  return_object <- list(init_SAD_par,init_curve_par,init_pro)
  names(return_object)<-c("init_SAD_par","init_curve_par","init_pro")
  return(return_object)
}
get_cluster <- function(data,k,input){
  Delta <- 10000; iter <- 1; itermax <- 100;
  
  get_biSAD1 <- function(par){
    n1=265;n2=146;n3=135
    
    AR1_get_matrix <- function(par,n){
      var_cov<-matrix(0,n,n)
      
      for (i in 1:n){
        for (j in 1:n){
          var_cov[i,j]<- par[2]^abs(j-i)*par[1]
        }
      }
      return(var_cov)
    }
    sig1 <- AR1_get_matrix(par[1:2],n1)
    sig2 <- AR1_get_matrix(par[3:4],n2)
    sig3 <- AR1_get_matrix(par[5:6],n3)
    
    sig12 <- array(0, dim=c(n1,n2))
    sig21 <- array(0, dim=c(n2,n1))
    sigma1 <- cbind(sig1,sig12)
    sigma2 <- cbind(sig21,sig2)
    sigma12 <- rbind(sigma1,sigma2)
    
    sig123 <- array(0, dim=c(n1+n2,n3))
    sig321 <- array(0, dim=c(n3,n1+n2))
    sigma1 <- cbind(sigma12,sig123)
    sigma2 <- cbind(sig321,sig3)
    sigma123 <- rbind(sigma1,sigma2)
    
    return(sigma123)
    
    return(sigma12)
  } 
  power_equation <- function(par,x){
    y <- par[1]*x^par[2]
    y
  }
  mle <- function(par,data,prob){
    par1 <- par[1:6]
    par2 <- matrix(par[-c(1:6)],nrow = k)
    miu <- t( sapply(1:k, function(c)c(power_equation(par2[c,1:2],meta_index1),
                                       power_equation(par2[c,3:4],meta_index2),
                                       power_equation(par2[c,5:6],meta_index3))))
    temp_S <- sapply(1:k,function(c)dmvnorm(data,
                                            miu[c,],
                                            get_biSAD1(par1))*prob[c])
    
    LL <- sum(-log(rowSums(temp_S)))
    
    return(LL)
  }
  
  cat(paste0("Start Clu Calculation ","\n","Cluster_number=",k))
  while ( Delta > 1e-3 && iter <= itermax ) {
    # initiation
    if(iter == 1){
      init_SAD_par <- input[[1]]
      init_curve_par <- input[[2]]
      pro <- input[[3]]
    }
    #E step, calculate the posterior probability
    old_par <- c(init_SAD_par,init_curve_par)
    LL_mem <- mle(old_par,data,pro)
    
    miu <- t( sapply(1:k, function(c)c(power_equation(init_curve_par[c,1:2],meta_index1),
                                       power_equation(init_curve_par[c,3:4],meta_index2),
                                       power_equation(init_curve_par[c,5:6],meta_index3))))
    mvn.c <- sapply(1:k, function(c) dmvnorm(data,
                                             miu[c,],
                                             get_biSAD1(init_SAD_par))*pro[c] )
    omega <- mvn.c/rowSums(mvn.c)
    
    #M step, calculate parameters
    pro <- colSums(omega)/sum(omega)
    
    new_par <- try(optim(old_par, mle, data=data, prob=pro, method = "Nelder-Mead"))
    if ('try-error' %in% class(new_par))
      break
    L_Value <- new_par$value
    init_SAD_par <- new_par$par[1:6]
    init_curve_par <- matrix(new_par$par[-c(1:6)],nrow = k)
    Delta <- abs(L_Value-LL_mem)
    
    cat("iter=",iter,"LL=",L_Value,'\n')
    iter <- iter+1; LL_mem <- L_Value
  } 
  
  BIC <- 2*(L_Value)+log(nrow(data))*length(old_par)
  return(BIC)
}



### net;subnet;age net
get_interaction <- function(data,col){

  n <- nrow(data)
  clean_data <- data
  gene_list <- list()
  m <- clean_data[,col]
  M <- clean_data[,-col]
  x_matrix <- M
  x_matrix <- as.matrix(x_matrix)
  name <- colnames(clean_data)
  
  ridge1_cv <- cv.glmnet(x = x_matrix, y = m,type.measure = "mse",family="gaussian",nfold = 10,alpha = 0)
  best_ridge_coef <- abs(as.numeric(coef(ridge1_cv, s = ridge1_cv$lambda.min))[-1])
  
  fit_res <- cv.glmnet(x = x_matrix, y = m,type.measure = "mse", family="gaussian",
                       nfold = 10,alpha = 1,
                       penalty.factor = 1/best_ridge_coef,
                       keep = TRUE)
  best_alasso_coef1 <- coef(fit_res, s = fit_res$lambda.min)
  
  return_obj = list(ind.name = name[col],
                    dep.name = best_alasso_coef1@Dimnames[[1]][best_alasso_coef1@i[-1]+1],
                    coefficient = best_alasso_coef1@x[-1])
  
  return(return_obj)
}
get_legendre_matrix <- function(x,legendre_order){
  legendre_coef <- legendre.polynomials(n = legendre_order, normalized=F)
  legendre_matrix <- as.matrix(as.data.frame(polynomial.values(
    polynomials = legendre_coef, x = scaleX(x, u = -1, v = 1))))
  colnames(legendre_matrix) <- paste0("legendre_",0:legendre_order)
  return(legendre_matrix[,2:(legendre_order+1)])
}
get_legendre_par <- function(y,legendre_order,x){
  #lm_method
  legendre_par <- as.numeric(coef(lm(y~get_legendre_matrix(x,legendre_order))))
  return(legendre_par)
}
legendre_fit <- function(par,x){
  legendre_order = length(par)
  fit <- sapply(1:length(par), function(c)
    par[c]*legendre.polynomials(n=legendre_order, normalized=F)[[c]])
  legendre_fit <- as.matrix(as.data.frame(polynomial.values(
    polynomials = fit, x = scaleX(x, u = -1, v = 1))))
  x_interpolation <- rowSums(legendre_fit)
  return(x_interpolation)
}
qdODEmod <- function(Time, State, Pars, power_par) {
  nn = length(Pars)
  ind_effect = paste0("alpha","*",names(State)[1])
  dep_effect = sapply(2:nn, function(c) paste0(paste0("beta",c-1),"*",names(State)[c]))
  dep_effect = paste0(dep_effect, collapse = "+")
  all_effect = paste0(ind_effect, "+", dep_effect)
  expr = parse(text = all_effect)
  
  with(as.list(c(State, Pars)), {
    dx = eval(expr)
    dy <- power_par[,1]*power_par[,2]*Time^(power_par[,2]-1)
    dind = alpha*x
    for(i in c(1:(nn-1))){
      tmp = paste0(paste0("beta",i),"*",paste0("y",i))
      expr2 = parse(text = tmp)
      assign(paste0("ddep",i),eval(expr2))
    }
    return(list(c(dx, dy, dind, mget(paste0("ddep",1:(nn-1))))))
  })
}
qdODE_ls <- function(pars, data, Time, power_par){
  n = length(pars)
  power_par = as.matrix(power_par)
  if (n==2) {
    Pars = c(alpha = pars[1], beta1 = pars[2:n])
    power_par = t(power_par)
    State = c(x=data[1,1],y1 = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  } else{
    Pars = c(alpha = pars[1], beta = pars[2:n])
    State = c(x=data[1,1],y = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  }
  out = as.data.frame(ode(func = qdODEmod, y = State, parms = Pars,
                          times = Time, power_par = power_par))
  X = as.numeric(data[1,])
  fit = as.numeric(out[,2])
  penalty = sum((out$ind[out$ind<0])^2)
  sse = crossprod(X-fit) + penalty 
  return(sse)
}
qdODE_fit <- function(pars, data, Time, power_par, LOP_order = 3, new_time = NULL, n_expand = 100){
  n = length(pars)
  if (n==2) {
    Pars = c(alpha = pars[1], beta1 = pars[2:n])
    power_par = t(power_par)
    State = c(x=data[1,1],y1 = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  } else{
    Pars = c(alpha = pars[1], beta = pars[2:n])
    State = c(x=data[1,1],y = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  }
  out = as.data.frame(ode(func = qdODEmod, y = State, parms = Pars,
                          times = Time, power_par = power_par))
  out2 = data.frame(x = out[,1], y = as.numeric(data[1,]), y.fit = out[,2],
                    ind = out[,(n+2)], dep = out[,(n+3):(ncol(out))])
  colnames(out2)[4:ncol(out2)] = c(rownames(data)[1], rownames(data)[2:n])
  rownames(out2) = NULL
  
  all_LOP_par = sapply(2:ncol(out2),function(c)get_legendre_par(out2[,c], LOP_order, out2$x))
  
  if (is.null(new_time)) {
    time2 = seq(min(Time), max(Time), length = n_expand)
    out3 = apply(all_LOP_par, 2, legendre_fit, x = time2)
    out3 = cbind(time2, out3)
  } else{
    out3 = apply(all_LOP_par, 2, legendre_fit, x = new_time)
    out3 = cbind(new_time, out3)
  }
  colnames(out3) = colnames(out2)
  result = list(fit = out2,
                predict = data.frame(out3),
                LOP_par = all_LOP_par)
  return(result)
}
qdODE_all <- function(result, relationship, Time = Time, power_par = power_par, i, init_pars = 1, LOP_order = 3, method = "ls",
                      new_time = NULL, n_expand = 100, maxit = 1e5){
  
  variable = c(relationship[[i]]$ind.name, relationship[[i]]$dep.name)
  data = result[variable,]
  
  if (length(variable)<=1) {
    qdODE.est = NA
    result = NA
    return.obj <- append(result, list(ODE.value = NA,
                                      parameters = NA))
  } else{
    power_par = power_par[variable,][-1,]
    n = nrow(data)
    pars_int = c(init_pars,relationship[[i]]$coefficient)
    if (method == "ls") {
      qdODE.est <- optim(pars_int, qdODE_ls, data = data, Time = Time, power_par = power_par,
                         method = "Nelder-Mead",
                         #method = "BFGS",
                         control = list(trace = TRUE, maxit = maxit))
      
      result <- qdODE_fit(pars = qdODE.est$par,
                          data = data,
                          power_par = power_par,
                          Time = Time)
      return.obj <- append(result, list(ODE.value = qdODE.est$value,
                                        parameters = qdODE.est$par))
    } else{
      qdODE.est <- optim(pars_int, qdODE_ls, data = data, Time = Time, power_par = power_par,
                         method = "Nelder-Mead",
                         control = list(trace = TRUE, maxit = maxit))
      
      result <- qdODE_fit(pars = qdODE.est$par,
                          data = data,
                          power_par = power_par,
                          Time = Time)
      return.obj <- append(result, list(ODE.value = qdODE.est$value,
                                        parameters = qdODE.est$par))
    }
  }
  return(return.obj)
}
qdODE_parallel <- function(result,relationship,Time,power_par, reduction = FALSE, thread = 2, maxit = 1e4){
  
  core.number <- thread
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterEvalQ(cl, {require(orthopolynom)})
  clusterEvalQ(cl, {require(deSolve)})
  clusterExport(cl, c("qdODEmod", "qdODE_ls", "qdODE_fit", "qdODE_all","get_legendre_matrix",
                      "get_legendre_par","legendre_fit","result","relationship","Time","power_par","maxit"), envir=environment())
  result = pblapply(1:nrow(result),function(c) qdODE_all(result = result,
                                                         relationship = relationship,
                                                         Time = Time,
                                                         power_par = power_par,
                                                         i = c,
                                                         maxit = maxit), cl = cl)
  stopCluster(cl)
  names(result) = rownames(par1)
  names(relationship) = rownames(par1)
  return_obj <- list(ode_result = result,
                     relationship = relationship)
  return(return_obj)
}

