
### Functions operating on single-pathogen Bayesian P-spline stan model fits

#' Function for estimating the model posterior of case incidence (median and credible intervals) for the single-pathogen Bayesian P-spline model.
#'
#' @param ps_fit the stan fit of the single-pathogen Bayesian P-spline model.
#' @param X list of dates/days over which the time series is defined.
#' @param num_days the number of dates/days (default is length of X).
#' @param time_labels descriptive labels for the values of X (e.g. may wish to provide the as.Date() objects to make plotting easier later).
#' @param days_per_knot the number of days between each knot (default of 5 days is used in all analysis).
#' @param spline_degree the spline degree used for the model (default of degree 3 is used in all analysis).
#'
#' @return data.frame() of modelled case incidence
#'
ps_single_incidence <- function(ps_fit,
                                X,
                                num_days = length(X), time_labels,
                                days_per_knot = 5, spline_degree = 3){
  
  
  X <- as.numeric(X)
  
  knots <- get_knots(X, days_per_knot = days_per_knot, spline_degree = spline_degree)
  num_knots <- length(knots)
  
  
  num_basis <- num_knots + spline_degree - 1
  num_data <- length(X)
  
  B_true <- splines::bs(seq(knots[1], knots[length(knots)], 1), knots = knots[2:(length(knots)-1)], degree=spline_degree, intercept = TRUE)
  B_true <- t(predict(B_true, X))
  
  
  post <- rstan::extract(ps_fit)
  
  a <- array(data=NA, dim=c(nrow(post$a), num_days))
  
  for(k in 1:nrow(post$a)){
    a[k,] <- as.array(post$a[k,]) %*% B_true
  }
  
  df <- data.frame()
  
  for(i in 1:num_days){
    total <- matrix(data=0, nrow=1, ncol=nrow(a))
    
    quan<- quantile(a[,i], c(0.5,0.025, 0.25, 0.75, 0.975))
    row <- data.frame(time = time_labels[i],
                      t_step = i,
                      y=exp(quan[[1]]),
                      lb_50 = exp(quan[[3]]),
                      ub_50 = exp(quan[[4]]),
                      lb_95 = exp(quan[[2]]),
                      ub_95 = exp(quan[[5]]))
    df <- rbind(df, row)
    
  }
  
  df
  
  
}


#' Function for estimating the model posterior of case incidence including day-of-the-week effects (median and credible intervals) for the single-pathogen Bayesian P-spline model.
#'
#' @param ps_fit the stan fit of the single-pathogen Bayesian P-spline model.
#' @param X list of dates/days over which the time series is defined.
#' @param DOW numerical values describing the group of each value in the X time series if distinct effects were included for distinct days (i.e. repeating values of 1:7 when each day of week has a seperate effect)
#' @param week_effect number descirbing the number of days that are modelled as having distinct effects (i.e. week_effect=7 when each day of week has a seperate effect)
#' @param num_days the number of dates/days (default is length of X).
#' @param time_labels descriptive labels for the values of X (e.g. may wish to provide the as.Date() objects to make plotting easier later).
#' @param days_per_knot the number of days between each knot (default of 5 days is used in all analysis).
#' @param spline_degree the spline degree used for the model (default of degree 3 is used in all analysis).
#'
#' @return data.frame() of modelled case incidence
#'
ps_single_incidence_dow <- function(ps_fit,
                                X, DOW, week_effect,
                                num_days = length(X), time_labels,
                                days_per_knot = 5, spline_degree = 3){
  
  
  X <- as.numeric(X)
  
  knots <- get_knots(X, days_per_knot = days_per_knot, spline_degree = spline_degree)
  num_knots <- length(knots)
  
  
  num_basis <- num_knots + spline_degree - 1
  num_data <- length(X)
  
  B_true <- splines::bs(seq(knots[1], knots[length(knots)], 1), knots = knots[2:(length(knots)-1)], degree=spline_degree, intercept = TRUE)
  B_true <- t(predict(B_true, X))
  
  
  post <- rstan::extract(ps_fit)
  
  a <- array(data=NA, dim=c(nrow(post$a), num_days))
  
  for(k in 1:nrow(post$a)){
    a[k,] <- as.array(post$a[k,]) %*% B_true
  }
  
  df <- data.frame()
  
  for(i in 1:num_days){
    total <- matrix(data=0, nrow=1, ncol=nrow(a))
    
    quan<- quantile(exp(a[,i])* week_effect*post$day_of_week_simplex[,DOW[i]], c(0.5,0.025, 0.25, 0.75, 0.975))
    row <- data.frame(time = time_labels[i],
                      t_step = i,
                      y=quan[[1]],
                      lb_50 = quan[[3]],
                      ub_50 = quan[[4]],
                      lb_95 = quan[[2]],
                      ub_95 = quan[[5]])
    df <- rbind(df, row)
    
  }
  
  df
  
  
}



#' Function for estimating the model posterior of growth rates (median and credible intervals) for the single-pathogen Bayesian P-spline model.
#'
#' @param ps_fit the stan fit of the single-pathogen Bayesian P-spline model.
#' @param X list of dates/days over which the time series is defined.
#' @param num_days the number of dates/days (default is length of X).
#' @param time_labels descriptive labels for the values of X (e.g. may wish to provide the as.Date() objects to make plotting easier later).
#' @param days_per_knot the number of days between each knot (default of 5 days is used in all analysis).
#' @param spline_degree the spline degree used for the model (default of degree 3 is used in all analysis).
#'
#' @return data.frame() of modelled growth rate
#'
ps_single_growth_rate <- function(ps_fit,
                           X,
                           num_days = length(X), time_labels,
                           days_per_knot = 5, spline_degree = 3){
  
  
  X <- as.numeric(X)
  
  knots <- get_knots(X, days_per_knot = days_per_knot, spline_degree = spline_degree)
  num_knots <- length(knots)
  
  
  num_basis <- num_knots + spline_degree - 1
  num_data <- length(X)
  
  B_true <- splines::bs(seq(knots[1], knots[length(knots)], 1), knots = knots[2:(length(knots)-1)], degree=spline_degree, intercept = TRUE)
  B_true <- t(predict(B_true, X))
  
  
  post <- rstan::extract(ps_fit)
  
  a <- array(data=NA, dim=c(nrow(post$a), num_days))
  
  for(k in 1:nrow(post$a)){
    a[k,] <- as.array(post$a[k,]) %*% B_true
  }
  
  df <- data.frame()
  
  for(i in 2:num_days){
    quan<- quantile(a[,i] - a[,(i-1)], c(0.5,0.025, 0.25, 0.75, 0.975))
    prop <- length(a[,i][ (a[,i] - a[,(i-1)])>0])/length(a[,i] - a[,(i-1)])
    row <- data.frame(time = time_labels[i],
                      t_step = i,
                      y=quan[[1]],
                      lb_50 = quan[[3]],
                      ub_50 = quan[[4]],
                      lb_95 = quan[[2]],
                      ub_95 = quan[[5]],
                      prop = prop)
    df <- rbind(df, row)
  
  }
  
  df

}



#' Function for estimating the model posterior of time-varying reproduction number (median and credible intervals) for the single-pathogen Bayesian P-spline model.
#'
#' @param ps_fit the stan fit of the single-pathogen Bayesian P-spline model.
#' @param X list of dates/days over which the time series is defined.
#' @param num_days the number of dates/days (default is length of X).
#' @param time_labels descriptive labels for the values of X (e.g. may wish to provide the as.Date() objects to make plotting easier later).
#' @param days_per_knot the number of days between each knot (default of 5 days is used in all analysis).
#' @param spline_degree the spline degree used for the model (default of degree 3 is used in all analysis).
#'
#' @return data.frame() of modelled time-varying reproduction number
#'
ps_single_Rt <- function(ps_fit,
                  X,
                  num_days = length(X), time_labels,
                  tau_max = 7,
                  gi_dist,
                  days_per_knot = 5, spline_degree = 3){
  
  
  X <- as.numeric(X)
  
  knots <- get_knots(X, days_per_knot = days_per_knot, spline_degree = spline_degree)
  num_knots <- length(knots)
  
  
  num_basis <- num_knots + spline_degree - 1
  num_data <- length(X)
  
  B_true <- splines::bs(seq(knots[1], knots[length(knots)], 1), knots = knots[2:(length(knots)-1)], degree=spline_degree, intercept = TRUE)
  B_true <- t(predict(B_true, X))
  
  
  post <- rstan::extract(ps_fit)
  
  a <- array(data=NA, dim=c(nrow(post$a), num_days))

  for(k in 1:nrow(post$a)){
    a[k,] <- as.array(post$a[k,]) %*% B_true
  }
  
  df <- data.frame()
  
  g_a <- sum(gi_dist(seq(0,tau_max-1,1)))
  
  for(i in tau_max:num_days){
      
      R_list <- matrix(0, nrow=length(a[,1]),ncol=1)
      
      for(k in 0:(tau_max-1)){
        R_list <- R_list + exp(a[,i-k])*gi_dist(k)
      }
      R_list <- exp(a[,i])/(R_list/g_a)
      
      quan<- quantile(R_list, c(0.5,0.025, 0.25, 0.75, 0.975))
      prop <- length(R_list[R_list>1]) / length(R_list)
      row <- data.frame(time = time_labels[i],
                        t_step = i,
                        y=quan[[1]],
                        lb_50 = quan[[3]],
                        ub_50 = quan[[4]],
                        lb_95 = quan[[2]],
                        ub_95 = quan[[5]],
                        prop = prop)
      df <- rbind(df, row)
    
  }
  
  df
  
  
  
}


