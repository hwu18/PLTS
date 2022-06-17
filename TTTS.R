args = commandArgs(trailingOnly=TRUE)

library(poisson)
library(MASS)
library(matrixcalc)
library(foreach)
library(doParallel)
likelihood = function(lambda, N_vec, infect_vec) {
  K = nrow(N_vec)
  num_infection = ncol(N_vec)
  indicator_mat = matrix(0, K, num_infection)
  for (i in 1:K){
    indicator_mat[i,] = as.numeric(infect_vec == i)
  }
  ans = 0
  temp = N_vec[1,]
  for (i in 1:(K-1)) {
    ans = ans + sum(indicator_mat[i+1,])*lambda[i]
    temp = temp + N_vec[i+1,]*exp(-lambda[i])
  }
  ans = ans + sum(log(temp))
  return(ans)
}

laplace_approx = function(N_vec, infect_vec) {
  K = nrow(N_vec)
  ans = optim(rep(0, K-1), likelihood, N_vec = N_vec, infect_vec = infect_vec, method = "BFGS",hessian = TRUE)
  pos_mode = ans$par
  cov_mat = tryCatch(solve(ans$hessian),
                     error = function(cond) {
                       return(matrix(0, K, K))
                     })
  use_MAP = tryCatch(min(eigen(cov_mat)$values) > 0 & sum(cov_mat) != 0,
                     error = function(cond) {
                       return(FALSE)
                     })
  use_MAP = !use_MAP
  return(list(mean = pos_mode, covariance = cov_mat, use_MAP = use_MAP)) 
}

draw_sample = function(cur_sol, n=1) {
  cur_mean = cur_sol$mean
  cur_cov = cur_sol$covariance
  use_MAP = cur_sol$use_MAP
  ans = if (use_MAP) t(replicate(n,cur_mean)) else mvrnorm(n, cur_mean, cur_cov)
  return(ans)
}

get_assignment_TS = function(cur_sol, num_new) {
  if (num_new == 0) {
    return(c())
  }
  data_mat = cbind(rep(0, num_new), draw_sample(cur_sol, num_new))
  ans = sapply(1:num_new, function(i) which.max(data_mat[i,]))
  return(ans)
}
TTTS = function(cur_sol) {
  I = which.max(c(0, draw_sample(cur_sol)))
  if (runif(1) < 0.5){
    assignment = I
  }else{
    J = which.max(c(0, draw_sample(cur_sol)))
    temp_count = 1
    while (J == I & temp_count < 30){
      temp_count = temp_count + 1
      J = which.max(c(0, draw_sample(cur_sol)))
    }
    assignment = J
  }
  return(assignment)
}

get_assignment_TTTS = function(cur_sol, n){
  return(replicate(n, TTTS(cur_sol)))
}

set.seed(1)
##baseline hazard
data = read.csv("case_daily_trends__united_states.csv", header = TRUE, skip = 3)
X = (0:349)
Y = rev(data$X7.Day.Moving.Avg[1:350])
fit = smooth.spline(X, Y)
h0 = function(t) {
  return((predict(fit, t)$y)/329484123)
}
plot_data = data.frame(x = 0:199, y = sapply(0:199, h0))
library(ggplot2)
p <- ggplot(plot_data, aes(x, 10000*y)) +
  geom_line(color="red", size=2)+theme(text = element_text(size=20)) + xlab("time") + ylab("number of infections per 10,000 people")

K = 6
T_max = 200
num_volunteer = 60000
hazard_vec = exp(c(0, -1.2, -1.5, -2.2, -2.4, -3))
T_vec = seq(from = 0, to = T_max, length.out = T_max+1)
B = 1000
num_new_volunteer = num_volunteer/T_max

run_TS = function(top_two) {
  count_stats_TS = matrix(0, K, T_max)
  infect_stats_TS = matrix(0, K, T_max)
  num_assign = matrix(0, K, B)
  for (j in 1:B) {
    ### TS
    print(j)
    num_count = numeric(K) #keep track of number of volunteers in each group
    total_count = numeric(K) #total number of volunteers in each group
    N_vec = matrix(0, K, 0)
    infect_vec = c()
    prev_time = 0
    for (q in 1:length(T_vec)) {
      time_step = T_vec[q]
      delta_t = time_step - prev_time
      cur_infection = sapply(1:K, function(i) rbinom(1, size = num_count[i], prob = delta_t*h0(prev_time)*hazard_vec[i]))
      if (q >= 2) {
        infect_stats_TS[,q-1] = infect_stats_TS[,q-1] + cur_infection
      }
      if (sum(cur_infection) > 0) {
        N_vec = cbind(N_vec, matrix(num_count, K, sum(cur_infection)))
        for (idx in 1:K) {
          infect_vec = c(infect_vec, rep(idx, cur_infection[idx]))
        }
        cur_sol = laplace_approx(N_vec, infect_vec)
      }
      
      if (length(infect_vec) == 0) {
        assign_vec = sample(1:K, num_new_volunteer, replace = TRUE, prob = rep(1,K)/K)
        assign_vec = sapply(1:K, function(i) sum(assign_vec==i))
        num_count = num_count + assign_vec
        total_count = total_count + assign_vec
        prev_time = time_step
        count_stats_TS[,time_step] = count_stats_TS[,time_step] + assign_vec
        next
      }
      num_count = num_count - cur_infection
      if (q < length(T_vec)) {
        if (top_two) {
          assign_vec = get_assignment_TTTS(cur_sol, num_new_volunteer)
        }else{
          assign_vec = get_assignment_TS(cur_sol, num_new_volunteer)
        }
        assign_vec = sapply(1:K, function(i) sum(assign_vec==i))
        num_count = num_count + assign_vec
        total_count = total_count + assign_vec
        prev_time = time_step
        count_stats_TS[,q] = count_stats_TS[,q] + assign_vec
      }
    }
    num_assign[,j] = total_count
  }
  return(list(count_stats_TS, infect_stats_TS, num_assign))
}

run_ans = run_TS(TRUE)
saveRDS(run_ans[[1]], "count_stats_TTTS.Rdata")
saveRDS(run_ans[[2]], "infect_stats_TTTS.Rdata")
saveRDS(run_ans[[3]], "count_each_run_TTTS.Rdata")
print("finished TTTS")