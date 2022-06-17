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

run_ans = run_TS(FALSE)
saveRDS(run_ans[[1]], "count_stats_TS.Rdata")
saveRDS(run_ans[[2]], "infect_stats_TS.Rdata")
saveRDS(run_ans[[3]], "count_each_run_TS.Rdata")
run_ans = run_TS(TRUE)
saveRDS(run_ans[[1]], "count_stats_TTTS.Rdata")
saveRDS(run_ans[[2]], "infect_stats_TTTS.Rdata")
saveRDS(run_ans[[3]], "count_each_run_TTTS.Rdata")
print("finished TTTS")

### DEW algorithm 
run_DEW = function(eta) {
  
  x = foreach (j=1:B) %dopar% {
    print(j)
    data_mat = matrix(0, num_volunteer, 4) #first column group index, second column enroll time, 
    #third column status (1 infected, 0 not), fourth column infected time 
    data_mat[,2] = rep(Inf, num_volunteer)
    data_mat[,3] = rep(-1, num_volunteer)
    data_mat[,4] = rep(Inf, num_volunteer)
    count_stats_DEW = matrix(0, K, T_max)
    infect_stats_DEW = matrix(0, K, T_max)
    
    weight_vec = rep(1, K)
    prob_mat = matrix(1, T_max, K)
    prev_time = 0
    for (q in 1:length(T_vec)) {
      time_step = T_vec[q]
      delta_t = time_step - prev_time
      if (q < length(T_vec)) {
        prob_mat[q,] = weight_vec/sum(weight_vec)
      }
      enrolled_idx = if (q == 1) c() else which(data_mat[,2] < time_step & data_mat[,3] != 1)
      infection_vec = if (length(enrolled_idx) == 0) c() else sapply(enrolled_idx, function(i) rbinom(1, size = 1, prob = delta_t*h0(prev_time)*hazard_vec[data_mat[i,1]]))
      if (length(infection_vec) > 0) {
        infected_idx = enrolled_idx[which(infection_vec == 1)]
        if (q >= 2) {
          infect_stats_DEW[,q-1] = infect_stats_DEW[,q-1] + sapply(1:K, function(i) sum(data_mat[infected_idx,1] == i)) 
        }
        data_mat[infected_idx, 3] = 1
        data_mat[infected_idx, 4] = time_step
        loss_vec = numeric(K)
        for (i in infected_idx) {
          group_idx = data_mat[i,1]
          enroll_time = data_mat[i,2]
          loss_vec[group_idx] = loss_vec[group_idx] + 1/prob_mat[which(T_vec == enroll_time),group_idx]
        }
        weight_vec = weight_vec*exp(-eta*loss_vec)
      }
      if (q < length(T_vec)) {
        assignment = sample(1:K, size = num_new_volunteer, prob = prob_mat[q,], replace = TRUE)
        count_stats_DEW[,q] = count_stats_DEW[,q] + sapply(1:K, function(i) sum(assignment == i))
        data_mat[((q-1)*num_new_volunteer+1):(q*num_new_volunteer),1] = assignment
        data_mat[((q-1)*num_new_volunteer+1):(q*num_new_volunteer),2] = time_step
        data_mat[((q-1)*num_new_volunteer+1):(q*num_new_volunteer),3] = 0
        prev_time = time_step 
      }
    }
    list(count_stats_DEW, infect_stats_DEW)
  }
  return(x)
}

ncore <- as.integer(args[1])

print(paste('running on', ncore, 'cores'))
registerDoParallel(ncore)

DEW_answers = list()
eta_vec = c(0.2, 0.1, 0.03, 0.01)
for (i in 1:length(eta_vec)) {
  DEW_answers[[i]] = run_DEW(eta_vec[i])
  print("finished!")
}
saveRDS(DEW_answers, "DEW.Rdata")

## Best vaccine identification
K = 6
T_max = 200
num_vec = c(60000)
hazard_vec = exp(c(0, -1.2, -1.5, -2.2, -2.4, -3))
T_vec = seq(from = 0, to = T_max, length.out = T_max+1)
B = 1000

best_arm_TS = function(num_volunteer) {
  num_new_volunteer = num_volunteer/T_max
  x = foreach (j=1:B, .combine = c) %dopar% {
    ### TS
    num_count = numeric(K) #keep track of number of volunteers in each group
    total_count = numeric(K) #total number of volunteers in each group
    N_vec = matrix(0, K, 0)
    infect_vec = c()
    prev_time = 0
    for (q in 1:length(T_vec)) {
      time_step = T_vec[q]
      delta_t = time_step - prev_time
      cur_infection = sapply(1:K, function(i) rbinom(1, size = num_count[i], prob = delta_t*h0(prev_time)*hazard_vec[i]))
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
        next
      }
      num_count = num_count - cur_infection
      if (q < length(T_vec)) {
        assign_vec = get_assignment_TTTS(cur_sol, num_new_volunteer)
        assign_vec = sapply(1:K, function(i) sum(assign_vec==i))
        num_count = num_count + assign_vec
        total_count = total_count + assign_vec
        prev_time = time_step
      }
    }
    which.max(c(0, cur_sol$mean))
  }
  return(x)
}

best_arm_result_TS = list()
for (i in 1:length(num_vec)) {
  best_arm_result_TS[[i]] = best_arm_TS(num_vec[i])
}
saveRDS(best_arm_result_TS, "best_arm_TS.Rdata")
print("finished TS best arm!")

best_arm_DEW = function(num_volunteer, eta) {
  num_new_volunteer = num_volunteer/T_max
  x = foreach (j=1:B, .combine = c) %dopar% {
    print(j)
    data_mat = matrix(0, num_volunteer, 4) #first column group index, second column enroll time, 
    #third column status (1 infected, 0 not), fourth column infected time 
    data_mat[,2] = rep(Inf, num_volunteer)
    data_mat[,3] = rep(-1, num_volunteer)
    data_mat[,4] = rep(Inf, num_volunteer)
  
    weight_vec = rep(1, K)
    prob_mat = matrix(1, T_max, K)
    prev_time = 0
    for (q in 1:length(T_vec)) {
      time_step = T_vec[q]
      delta_t = time_step - prev_time
      if (q < length(T_vec)) {
        prob_mat[q,] = weight_vec/sum(weight_vec)
      }
      enrolled_idx = if (q == 1) c() else which(data_mat[,2] < time_step & data_mat[,3] != 1)
      infection_vec = if (length(enrolled_idx) == 0) c() else sapply(enrolled_idx, function(i) rbinom(1, size = 1, prob = delta_t*h0(prev_time)*hazard_vec[data_mat[i,1]]))
      if (length(infection_vec) > 0) {
        infected_idx = enrolled_idx[which(infection_vec == 1)]
        data_mat[infected_idx, 3] = 1
        data_mat[infected_idx, 4] = time_step
        loss_vec = numeric(K)
        for (i in infected_idx) {
          group_idx = data_mat[i,1]
          enroll_time = data_mat[i,2]
          loss_vec[group_idx] = loss_vec[group_idx] + 1/prob_mat[which(T_vec == enroll_time),group_idx]
        }
        weight_vec = weight_vec*exp(-eta*loss_vec)
      }
      if (q < length(T_vec)) {
        assignment = sample(1:K, size = num_new_volunteer, prob = prob_mat[q,], replace = TRUE)
        data_mat[((q-1)*num_new_volunteer+1):(q*num_new_volunteer),1] = assignment
        data_mat[((q-1)*num_new_volunteer+1):(q*num_new_volunteer),2] = time_step
        data_mat[((q-1)*num_new_volunteer+1):(q*num_new_volunteer),3] = 0
        prev_time = time_step 
      }
    }
    which.max(weight_vec)
  }
  return(x)
}

best_arm_result_DEW = list()
eta_vec = c(0.1, 0.01)
for (i in 1:length(num_vec)) {
  best_arm_result_DEW[[i]] = list()
  for (j in 1:length(eta_vec)) {
   best_arm_result_DEW[[i]][[j]] = best_arm_DEW(num_vec[i], eta_vec[j]) 
  }
}
saveRDS(best_arm_result_DEW, "best_arm_DEW.Rdata")
print("finished DEW best arm!")

best_arm_RCT = function(num_volunteer) {
  num_new_volunteer = num_volunteer/T_max
  x = foreach(j = 1:B, .combine = c) %dopar% {
    num_count = numeric(K)
    total_count = numeric(K)
    infect_vec = numeric(K)
    prev_time= 0
    
    for (time_step in T_vec) {
      delta_t = time_step - prev_time
      cur_infection = sapply(1:K, function(i) rbinom(1, size = num_count[i], prob = delta_t*h0(prev_time)*hazard_vec[i]))
      infect_vec = infect_vec + cur_infection
      num_count = num_count - cur_infection
      if (time_step == T_max + 1) {
        break
      }
      assign_vec = sample(1:K, num_new_volunteer, replace = TRUE, prob = rep(1,K)/K)
      assign_vec = sapply(1:K, function(i) sum(assign_vec==i))
      num_count = num_count + assign_vec
      total_count = total_count + assign_vec
      prev_time = time_step
    }
    min_idx = which(infect_vec/total_count == min(infect_vec/total_count))
    best_idx = if (length(min_idx) > 1) sample(min_idx, 1, prob = rep(1/length(min_idx), length(min_idx)))  else min_idx
    best_idx
  }
  return(x)
}

best_arm_result_RCT = list()
for (i in 1:length(num_vec)) {
  best_arm_result_RCT[[i]] = best_arm_RCT(num_vec[i])
}
saveRDS(best_arm_result_RCT, "best_arm_RCT.Rdata")
print("finished RCT best arm!")

stopImplicitCluster()
