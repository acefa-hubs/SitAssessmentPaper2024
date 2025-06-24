
### Functions operating on multi-pathogen Bayesian P-spline stan model fits

#' Function for getting the location of knots (required for fitting P-spline models)
#'
#' @param X list of dates/days over which the time series (for model fitting) is defined.
#' @param days_per_knot the number of days between each knot (default of 5 days is used in all analysis)
#' @param spline_degree the spline degree used for the model (default of degree 3 is used in all analysis)
#'
#' @return location of knots
#'
get_knots <- function(X, days_per_knot=5, spline_degree=3){
   
  X <- as.numeric(X)
  
  num_knots <- ceiling((max(X)-min(X))/days_per_knot)
  
  first_knot <- min(X) - spline_degree*days_per_knot
  final_knot <- first_knot + days_per_knot*num_knots + 2*spline_degree*days_per_knot
  
  knots <- seq(first_knot, final_knot, by=days_per_knot)
  
}


#' Function for estimating the model posterior of case incidence (median and credible intervals) for the multi-pathogen Bayesian P-spline model.
#'
#' @param ps_fit the stan fit of the multi-pathogen Bayesian P-spline model.
#' @param X list of dates/days over which the time series is defined.
#' @param num_days the number of dates/days (default is length of X).
#' @param time_labels descriptive labels for the values of X (e.g. may wish to provide the as.Date() objects to make plotting easier later).
#' @param days_per_knot the number of days between each knot (default of 5 days is used in all analysis.
#' @param spline_degree the spline degree used for the model (default of degree 3 is used in all analysis).
#' @param num_path the number of pathogens that were included in the multi pathogen Bayesian P-spline model (defaults to standard options for influenza-like illness).
#' @param pathogen_names the names of the pathogens that were included in the multi pathogen Bayesian P-spline model fit (defaults to standard options for influenza-like illness).
#'
#' @return data.frame() of modelled case incidence for each pathogen
#'
ps_incidence <- function(ps_fit,
                         X,
                         num_days = length(X), time_labels,
                         days_per_knot = 5, spline_degree = 3,
                         num_path=4, pathogen_names=c("Influenza A H3N2", "Influenza A H1N1", "Influenza B", "Unknown")){
  
  
  X <- as.numeric(X)
  
  knots <- get_knots(X, days_per_knot = days_per_knot, spline_degree = spline_degree)
  num_knots <- length(knots)
  
  
  num_basis <- num_knots + spline_degree - 1
  num_data <- length(X)
  
  B_true <- splines::bs(seq(knots[1], knots[length(knots)], 1), knots = knots[2:(length(knots)-1)], degree=spline_degree, intercept = TRUE)
  B_true <- t(predict(B_true, X))
  
  
  post <- rstan::extract(ps_fit)
  
  a <- array(data=NA, dim=c(nrow(post$a), num_path, num_days))
  
  for(j in 1:num_path){
    for(k in 1:nrow(post$a)){
      a[k,j,] <- as.array(post$a[k,j,]) %*% B_true
    }
  }
  
  df <- data.frame()
  
  for(i in 1:num_days){
    
    total <- matrix(data=0, nrow=1, ncol=nrow(a))
    
    for(j in 1:num_path){
      quan<- quantile(a[,j,i], c(0.5,0.025, 0.25, 0.75, 0.975))
      row <- data.frame(time = time_labels[i],
                        t_step = i,
                        y=exp(quan[[1]]),
                        lb_50 = exp(quan[[3]]),
                        ub_50 = exp(quan[[4]]),
                        lb_95 = exp(quan[[2]]),
                        ub_95 = exp(quan[[5]]),
                        pathogen = pathogen_names[j])
      df <- rbind(df, row)
      
      total <- total + exp(a[,j,i])
      
    }
    
    quan<- quantile(total, c(0.5,0.025, 0.25, 0.75, 0.975))
    row <- data.frame(time = time_labels[i],
                      t_step = i,
                      y=quan[[1]],
                      lb_50 = quan[[3]],
                      ub_50 = quan[[4]],
                      lb_95 = quan[[2]],
                      ub_95 = quan[[5]],
                      pathogen = "Total")
    df <- rbind(df, row)
    
  }
  
  df
  
  
}


#' Function for estimating the model posterior of case incidence including day-of-the-week effects (median and credible intervals) for the multi-pathogen Bayesian P-spline model.
#'
#' @param ps_fit the stan fit of the multi-pathogen Bayesian P-spline model.
#' @param X list of dates/days over which the time series is defined.
#' @param num_days the number of dates/days (default is length of X).
#' @param time_labels descriptive labels for the values of X (e.g. may wish to provide the as.Date() objects to make plotting easier later).
#' @param DOW numerical values describing the group of each value in the X time series if distinct effects were included for distinct days (i.e. repeating values of 1:7 when each day of week has a seperate effect)
#' @param week_effect number descirbing the number of days that are modelled as having distinct effects (i.e. week_effect=7 when each day of week has a seperate effect)
#' @param days_per_knot the number of days between each knot (default of 5 days is used in all analysis.
#' @param spline_degree the spline degree used for the model (default of degree 3 is used in all analysis).
#' @param num_path the number of pathogens that were included in the multi pathogen Bayesian P-spline model (defaults to standard options for influenza-like illness).
#' @param pathogen_names the names of the pathogens that were included in the multi pathogen Bayesian P-spline model fit (defaults to standard options for influenza-like illness).
#'
#' @return data.frame() of modelled case incidence including day-of-the-week effects for each pathogen
#'
ps_incidence_dow <- function(ps_fit,
                         X,
                         num_days = length(X), time_labels, DOW, week_effect,
                         days_per_knot = 5, spline_degree = 3,
                         num_path=4, pathogen_names=c("Influenza A H3N2", "Influenza A H1N1", "Influenza B", "Unknown")){
  
  
  X <- as.numeric(X)
  
  knots <- get_knots(X, days_per_knot = days_per_knot, spline_degree = spline_degree)
  num_knots <- length(knots)
  
  
  num_basis <- num_knots + spline_degree - 1
  num_data <- length(X)
  
  B_true <- splines::bs(seq(knots[1], knots[length(knots)], 1), knots = knots[2:(length(knots)-1)], degree=spline_degree, intercept = TRUE)
  B_true <- t(predict(B_true, X))
  
  
  post <- rstan::extract(ps_fit)
  
  a <- array(data=NA, dim=c(nrow(post$a), num_path, num_days))
  
  for(j in 1:num_path){
    for(k in 1:nrow(post$a)){
      a[k,j,] <- as.array(post$a[k,j,]) %*% B_true
    }
  }
  
  df <- data.frame()
  
  for(i in 1:num_days){
    
    total <- matrix(data=0, nrow=1, ncol=nrow(a))
    
    for(j in 1:num_path){
      quan<- quantile(exp(a[,j,i])* week_effect*post$day_of_week_simplex[,DOW[i]], c(0.5,0.025, 0.25, 0.75, 0.975))
      row <- data.frame(time = time_labels[i],
                        t_step = i,
                        y=quan[[1]],
                        lb_50 = quan[[3]],
                        ub_50 = quan[[4]],
                        lb_95 = quan[[2]],
                        ub_95 = quan[[5]],
                        pathogen = pathogen_names[j])
      df <- rbind(df, row)
      
      total <- total + exp(a[,j,i])
      
    }
    
    quan<- quantile(total*week_effect*post$day_of_week_simplex[,DOW[i]], c(0.5,0.025, 0.25, 0.75, 0.975))
    row <- data.frame(time = time_labels[i],
                      t_step = i,
                      y=quan[[1]],
                      lb_50 = quan[[3]],
                      ub_50 = quan[[4]],
                      lb_95 = quan[[2]],
                      ub_95 = quan[[5]],
                      pathogen = "Total")
    df <- rbind(df, row)
    
  }
  
  df
  
  
}





#' Function for estimating the model posterior of relative proportions of different pathogens (median and credible intervals) for the multi-pathogen Bayesian P-spline model.
#'
#' @param ps_fit the stan fit of the multi-pathogen Bayesian P-spline model.
#' @param X list of dates/days over which the time series is defined.
#' @param num_days the number of dates/days (default is length of X).
#' @param time_labels descriptive labels for the values of X (e.g. may wish to provide the as.Date() objects to make plotting easier later).
#' @param days_per_knot the number of days between each knot (default of 5 days is used in all analysis.
#' @param spline_degree the spline degree used for the model (default of degree 3 is used in all analysis).
#' @param num_path the number of pathogens that were included in the multi pathogen Bayesian P-spline model (defaults to standard options for influenza-like illness).
#' @param comb_num list of lists describing the pathogens to be included in the numerator of the proportions (deaults to standard options for influenza-like illness, see below)
#' @param comb_den list of lists describing the pathogens to be included in the denominator of the proportions (deaults to standard options for influenza-like illness, see below)
#' @param comb_names descriptive names descibing the different modelled relative proportions being estimated (deaults to standard options for influenza-like illness: prop influenza positive; prop influenza A; prop influenza B; prop H3N2 vs H1N1; prop H1N1 vs H3N2)
#'
#' @return data.frame() of modelled relative proportions of different pathogens (for comparison to data)
#'
ps_proportion <- function(ps_fit,
                         X,
                         num_days = length(X), time_labels,
                         days_per_knot = 5, spline_degree = 3,
                         num_path=4, 
                         comb_num=list(c(1,2,3), c(1,2), c(3), c(1), c(2)),
                         comb_den=list(c(1,2,3,4), c(1,2,3,4), c(1,2,3,4), c(1,2), c(1,2) ),
                         comb_names=c("Influenza", "Influenza A", "Influenza B", "H3N2", "H1N1") ){
  
  
  X <- as.numeric(X)
  
  knots <- get_knots(X, days_per_knot = days_per_knot, spline_degree = spline_degree)
  num_knots <- length(knots)
  
  
  num_basis <- num_knots + spline_degree - 1
  num_data <- length(X)
  
  B_true <- splines::bs(seq(knots[1], knots[length(knots)], 1), knots = knots[2:(length(knots)-1)], degree=spline_degree, intercept = TRUE)
  B_true <- t(predict(B_true, X))
  
  
  post <- rstan::extract(ps_fit)
  
  a <- array(data=NA, dim=c(nrow(post$a), num_path, num_days))
  
  for(j in 1:num_path){
    for(k in 1:nrow(post$a)){
      a[k,j,] <- as.array(post$a[k,j,]) %*% B_true
    }
  }
  
  
  
  
  df <- data.frame()
  
  for(i in 1:num_days){
    
    total <- matrix(data=0, nrow=1, ncol=nrow(a))
    
    for(j in 1:length(comb_num)){
      num_index <- comb_num[[j]]
      den_index <- comb_den[[j]]
      
      if(length(num_index)>1){
        num <- rowSums(exp(a[,num_index,i]))
      } else{
        num <- exp(a[,num_index,i])
      }
      
      den <- rowSums(exp(a[,den_index,i]))
      
      quan<- quantile(num/den, c(0.5,0.025, 0.25, 0.75, 0.975))
      row <- data.frame(time = time_labels[i],
                        t_step = i,
                        y=quan[[1]],
                        lb_50 = quan[[3]],
                        ub_50 = quan[[4]],
                        lb_95 = quan[[2]],
                        ub_95 = quan[[5]],
                        pathogen = comb_names[j])
      df <- rbind(df, row)
      
    }
    
  }
  
  df
  
  
}



#' Function for estimating the model posterior of growth rates (median and credible intervals) for the multi-pathogen Bayesian P-spline model.
#'
#' @param ps_fit the stan fit of the multi-pathogen Bayesian P-spline model.
#' @param X list of dates/days over which the time series is defined.
#' @param num_days the number of dates/days (default is length of X).
#' @param time_labels descriptive labels for the values of X (e.g. may wish to provide the as.Date() objects to make plotting easier later).
#' @param days_per_knot the number of days between each knot (default of 5 days is used in all analysis.
#' @param spline_degree the spline degree used for the model (default of degree 3 is used in all analysis).
#' @param num_path the number of pathogens that were included in the multi pathogen Bayesian P-spline model (defaults to standard options for influenza-like illness).
#' @param pathogen_names the names of the pathogens that were included in the multi pathogen Bayesian P-spline model fit (defaults to standard options for influenza-like illness).
#'
#' @return data.frame() of modelled growth rates for each pathogen
#'
ps_growth_rate <- function(ps_fit,
                           X,
                           num_days = length(X), time_labels,
                           days_per_knot = 5, spline_degree = 3,
                           num_path=4, pathogen_names=c("Influenza A H3N2", "Influenza A H1N1", "Influenza B", "Unknown")){
  
  
  X <- as.numeric(X)
  
  knots <- get_knots(X, days_per_knot = days_per_knot, spline_degree = spline_degree)
  num_knots <- length(knots)
  
  
  num_basis <- num_knots + spline_degree - 1
  num_data <- length(X)
  
  B_true <- splines::bs(seq(knots[1], knots[length(knots)], 1), knots = knots[2:(length(knots)-1)], degree=spline_degree, intercept = TRUE)
  B_true <- t(predict(B_true, X))
  
  
  post <- rstan::extract(ps_fit)
  
  a <- array(data=NA, dim=c(nrow(post$a), num_path, num_days))
  
  for(j in 1:num_path){
    for(k in 1:nrow(post$a)){
      a[k,j,] <- as.array(post$a[k,j,]) %*% B_true
    }
  }
  
  df <- data.frame()
  
  for(i in 2:num_days){
    
    total_i   <- matrix(data=0, nrow=1, ncol=nrow(a))
    total_im1 <- matrix(data=0, nrow=1, ncol=nrow(a))
    
    for(j in 1:num_path){
      quan<- quantile(a[,j,i] - a[,j,(i-1)], c(0.5,0.025, 0.25, 0.75, 0.975))
      prop <- length(a[,j,i][ (a[,j,i] - a[,j,(i-1)])>0])/length(a[,j,i] - a[,j,(i-1)])
      row <- data.frame(time = time_labels[i],
                        t_step = i,
                        y=quan[[1]],
                        lb_50 = quan[[3]],
                        ub_50 = quan[[4]],
                        lb_95 = quan[[2]],
                        ub_95 = quan[[5]],
                        pathogen = pathogen_names[j],
                        prop = prop)
      df <- rbind(df, row)
      
      total_i   <- total_i   + exp(a[,j,i])
      total_im1 <- total_im1 + exp(a[,j,(i-1)])
      
    }
    
    quan<- quantile(log(total_i)-log(total_im1), c(0.5,0.025, 0.25, 0.75, 0.975))
    prop <- length(total_i[ (log(total_i)-log(total_im1)) >0])/length(log(total_i)-log(total_im1))
    row <- data.frame(time = time_labels[i],
                      t_step = i,
                      y=quan[[1]],
                      lb_50 = quan[[3]],
                      ub_50 = quan[[4]],
                      lb_95 = quan[[2]],
                      ub_95 = quan[[5]],
                      pathogen = "Total",
                      prop = prop)
    df <- rbind(df, row)
    
  }
  
  df
  
  
  
}


#' Function for estimating the model posterior of time-varying reproduction numbers (median and credible intervals) for the multi-pathogen Bayesian P-spline model.
#'
#' @param ps_fit the stan fit of the multi-pathogen Bayesian P-spline model.
#' @param X list of dates/days over which the time series is defined.
#' @param num_days the number of dates/days (default is length of X).
#' @param time_labels descriptive labels for the values of X (e.g. may wish to provide the as.Date() objects to make plotting easier later).
#' @param tau_max the maximum number of previous days to consider when calculating Rt from the modelled case time-series: R(t) = I(t) / integral_0^tau_max( I(t-tau) g(tau) dtau )
#' @param gi_dist function that returns the value of the generation interval for different value of tau (the ammount of time before time t - see above equation)
#' @param days_per_knot the number of days between each knot (default of 5 days is used in all analysis.
#' @param spline_degree the spline degree used for the model (default of degree 3 is used in all analysis).
#' @param num_path the number of pathogens that were included in the multi pathogen Bayesian P-spline model (defaults to standard options for influenza-like illness).
#' @param pathogen_names the names of the pathogens that were included in the multi pathogen Bayesian P-spline model fit (defaults to standard options for influenza-like illness).
#'
#' @return data.frame() of modelled time-varying reproduction number for each pathogen
#'
ps_Rt <- function(ps_fit,
                  X,
                  num_days = length(X), time_labels,
                  tau_max = 7,
                  gi_dist,
                  days_per_knot = 5, spline_degree = 3,
                  num_path=4, pathogen_names=c("Influenza A H3N2", "Influenza A H1N1", "Influenza B", "Unknown")){
  
  
  X <- as.numeric(X)
  
  knots <- get_knots(X, days_per_knot = days_per_knot, spline_degree = spline_degree)
  num_knots <- length(knots)
  
  
  num_basis <- num_knots + spline_degree - 1
  num_data <- length(X)
  
  B_true <- splines::bs(seq(knots[1], knots[length(knots)], 1), knots = knots[2:(length(knots)-1)], degree=spline_degree, intercept = TRUE)
  B_true <- t(predict(B_true, X))
  
  
  post <- rstan::extract(ps_fit)
  
  a <- array(data=NA, dim=c(nrow(post$a), num_path, num_days))
  
  for(j in 1:num_path){
    for(k in 1:nrow(post$a)){
      a[k,j,] <- as.array(post$a[k,j,]) %*% B_true
    }
  }
  
  df <- data.frame()
  
  g_a <- sum(gi_dist(seq(0,tau_max-1,1)))
  
  for(i in tau_max:num_days){
    
    R_list_T <- matrix(0, nrow=length(a[,1,1]),ncol=1) #
    for(j in 1:num_path){
      
      R_list <- matrix(0, nrow=length(a[,1,1]),ncol=1)
      
      for(k in 0:(tau_max-1)){
        R_list <- R_list + exp(a[,j,i-k])*gi_dist(k)
        
        R_list_T <- R_list_T + exp(a[,j,i-k])*gi_dist(k)#
      }
      R_list <- exp(a[,j,i])/(R_list/g_a)
      
      quan<- quantile(R_list, c(0.5,0.025, 0.25, 0.75, 0.975))
      prop <- length(R_list[R_list>1]) / length(R_list)
      row <- data.frame(time = time_labels[i],
                        t_step = i,
                        y=quan[[1]],
                        lb_50 = quan[[3]],
                        ub_50 = quan[[4]],
                        lb_95 = quan[[2]],
                        ub_95 = quan[[5]],
                        pathogen = pathogen_names[j],
                        prop = prop)
      df <- rbind(df, row)
      
    }
    
    R_list_T <- rowSums(exp(a[ , ,i])) /(R_list_T/g_a)
    quan<- quantile(R_list_T, c(0.5,0.025, 0.25, 0.75, 0.975))
    prop <- length(R_list_T[R_list_T>1]) / length(R_list_T)
    row <- data.frame(time = time_labels[i],
                      t_step = i,
                      y=quan[[1]],
                      lb_50 = quan[[3]],
                      ub_50 = quan[[4]],
                      lb_95 = quan[[2]],
                      ub_95 = quan[[5]],
                      pathogen = "Total",
                      prop = prop)
    df <- rbind(df, row)
    
  }
  
  df
  
  
  
}


