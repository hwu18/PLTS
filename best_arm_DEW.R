args = commandArgs(trailingOnly=TRUE)

library(poisson)
library(MASS)
library(matrixcalc)
library(foreach)
library(doParallel)

set.seed(1)
##baseline hazard
data = read.csv("case_daily_trends__united_states.csv", header = TRUE, skip = 3)
X = (0:349)
Y = rev(data$X7.Day.Moving.Avg[1:350])
fit = smooth.spline(X, Y)
h0 = function(t) {
  return((predict(fit, t)$y)/329484123)
}
K = 6
T_max = 200
num_volunteer = 60000
hazard_vec = exp(c(0, -1.2, -1.5, -2.2, -2.4, -3))
T_vec = seq(from = 0, to = T_max, length.out = T_max+1)
B = 1000
num_new_volunteer = num_volunteer/T_max

ncore <- as.integer(args[1])

print(paste('running on', ncore, 'cores'))
registerDoParallel(ncore)

## Best vaccine identification
K = 6
T_max = 200
num_vec = c(60000)
hazard_vec = exp(c(0, -1.2, -1.5, -2.2, -2.4, -3))
T_vec = seq(from = 0, to = T_max, length.out = T_max+1)
B = 1000

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
eta_vec = c(1.6, 0.8, 0.4, 0.2, 0.1, 0.03, 0.01)
for (j in 1:length(eta_vec)) {
  best_arm_result_DEW[[j]] = best_arm_DEW(60000, eta_vec[j]) 
}
saveRDS(best_arm_result_DEW, "best_arm_DEW_60000.Rdata")
print("finished DEW best arm!")
stopImplicitCluster()