#=====================#
# DISCLAIMER
# This  code is provided "as is" and "with all faults" for educational and research use only.  
#   The author makes no representations or warranties of any kind concerning its correctness or 
#   suitability for any particular purposes.  Some numbers and graphs in published papers were
#   produced in Excel, not this code. If you find errors or ways to improve the code, please 
#   let me know at  R.G.ThomasNOSPAM AT kent.ac.uk
#=====================#

# Setting parameters 

R <- 2 # number of times drawdown traverse from 1 to b (top to bottom)

S <- 1.0 #stock price at time t
b <- 0.5 #barrier - for b = K case, set b just slightly less than K, to avoid divisions by zero problems.
tau <- 10 #time to maturity T - t (in years) 
r <- 0.015 #risk-free annual interest rate (convenient to set to zero)
q <- 0.01 #deferment rate (yield) (needs slightly different from r, to avoid division-by-zero problems in theta)
sigma <- 0.30 #annual volatility of the stock price (standard deviation)

set.seed(1930) #set the seed (if not set, it's taken from the computer clock)
N <- 6300 #N is number of time steps, e.g. 252 x 25 = 6300 steps for every working day over 25 years. 
nSim <- 100 #number of simulations (paths) 



dt <- tau/N #length of each time sub interval
Z <-  matrix(rnorm(nSim*N, mean=0, sd=1),nrow = nSim, ncol = N) #standard normal sample of N elements
dW <- Z*sqrt(dt) #Brownian motion increments (nSim simulations) x (N increments)
X <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
b_trailing <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
B <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
B_running_max <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Y_running_max <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Y <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Y_T <- vector(length=nSim) 

b_as_matrix <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
B_old <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
B_running_max_old <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Y_old <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Y_running_max_old <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)

Z_t_over_M_t <-matrix(numeric(nSim*N), nrow = nSim, ncol = N)

B_old_double_top <- array(dim=c(nSim,N,R)) # extra s-dimension is to give loops for recursive reflection
B_old_double_bottom <- array(dim=c(nSim,N,R))
B_running_min_old_double <- array(dim=c(nSim,N,R))
B_running_max_old_double <- array(dim=c(nSim,N,R))
Y_old_drawdown <- array(dim=c(nSim,N,R))
Y_running_min_old_double <- array(dim=c(nSim,N,R))
Y_running_max_old_double <- array(dim=c(nSim,N,R))

X[,1] <- S
Y[,1] <- S
b_trailing[,1] <- b
B_running_max[,1] <- b_trailing[,1]/X[,1]
Y_running_max[,1] <- S
Z_t_over_M_t[,1] <- 1

Y_old[,1] <- X[,1]
b_as_matrix[,] <- b
B_running_max_old[,1] <- b_as_matrix[,1]/X[,1]
Y_running_max_old[,1] <- S

Y_old_drawdown[,1,] <- X[,1]
B_running_min_old_double[,1,] <- 1 /X[,1]
B_running_max_old_double[,1,] <- b_as_matrix[,1]/X[,1]
Y_running_min_old_double[,1,] <- S
Y_running_max_old_double[,1,] <- S

#Then do the loop: X is GBM, Y is reflection at trailing barrier, Y_old is reflection at fixed barrier

for(j in 2:N){
    
    X[,j] <- X[,j-1]*exp((r - q -0.5*sigma^2)*dt + sigma*dW[,j])
    B[,j] <- b_trailing[,j-1]/X[,j] #Have to use b_trailing[j-1] because don't know new b_trailing before calculating latest Y
    B_running_max[,j] <- ifelse(B[,j] > B_running_max[,j-1], B[,j], B_running_max[,j-1])
    Y[,j] <- X[,j]*pmax(1,B_running_max[,j])
    Y_running_max[,j] <- ifelse(Y[,j] > Y_running_max[,j-1], Y[,j], Y_running_max[,j-1])
    b_trailing[,j] <- ifelse(Y[,j] > Y_running_max[,j-1], b * Y[,j], b_trailing[,j-1])
    Z_t_over_M_t[,j] <- Y[,j] / Y_running_max[,j]
      
    B_old[,j] <- b_as_matrix[,j]/X[,j] 
    B_running_max_old[,j] <- ifelse(B_old[,j] > B_running_max_old[,j-1], B_old[,j], B_running_max_old[,j-1])
    Y_old[,j] <- X[,j]*pmax(1,B_running_max_old[,j])
    Y_running_max_old[,j] <- ifelse(Y_old[,j] > Y_running_max_old[,j-1], Y_old[,j], Y_running_max_old[,j-1])
    

  }  

# Now construct the drawdown of the "GBM reflected at trailing barrier", by recursive reflection 
# of the original GBM at 1 and b.


for (j in 2:N){ # Reflect original GBM, X, at the top first:
  
B_old_double_top[,j,1] <- 1 /X[,j]  
B_running_min_old_double[,j,1] <- ifelse(B_old_double_top[,j,1] < B_running_min_old_double[,j-1,1], 
                                         B_old_double_top[,j,1], B_running_min_old_double[,j-1,1])
Y_old_drawdown[,j,1] <- X[,j] * pmin(1,B_running_min_old_double[,j,1])   
Y_running_min_old_double[,j,1] <- ifelse(Y_old_drawdown[,j,1] < Y_running_min_old_double[,j-1,1], 
                                         Y_old_drawdown[,j,1], Y_running_min_old_double[,j-1,1])

}

for (j in 2:N){ # Then reflect result at the bottom:

  B_old_double_bottom[,j,1] <- b_as_matrix[,j]/Y_old_drawdown[,j,1] # from above
  B_running_max_old_double[,j,1] <- ifelse(B_old_double_bottom[,j,1] > B_running_max_old_double[,j-1,1],
                                           B_old_double_bottom[,j,1], B_running_max_old_double[,j-1,1])
  Y_old_drawdown[,j,1] <- Y_old_drawdown[,j,1] * pmax(1,B_running_max_old_double[,j,1])
  Y_running_max_old_double[,j,1] <- ifelse(Y_old_drawdown[,j,1] > Y_running_max_old_double[,j-1,1],
                                           Y_old_drawdown[,j,1], Y_running_max_old_double[,j-1,1])

}

# Now reflect the result at {top, bottom} for another R-1 times
# (R is observed number of times the drawdown traverses 1 to b):


for (s in 2:R){

  for (j in 2:N){  # top reflection

    B_old_double_top[,j,s] <- 1 /Y_old_drawdown[,j,s-1]  # has to be s-1 here, not found for s yet
    B_running_min_old_double[,j,s] <- ifelse(B_old_double_top[,j,s] < B_running_min_old_double[,j-1,s],
                                             B_old_double_top[,j,s], B_running_min_old_double[,j-1,s])
    Y_old_drawdown[,j,s] <- Y_old_drawdown[,j,s-1] * pmin(1,B_running_min_old_double[,j,s])
    Y_running_min_old_double[,j,s] <- ifelse(Y_old_drawdown[,j,s] < Y_running_min_old_double[,j-1,s],
                                             Y_old_drawdown[,j,s], Y_running_min_old_double[,j-1,s])

  }

  for (j in 2:N){ # bottom reflection

    B_old_double_bottom[,j,s] <- b_as_matrix[,j]/Y_old_drawdown[,j,s] # s as final subscript here, using output from directly above
    B_running_max_old_double[,j,s] <- ifelse(B_old_double_bottom[,j,s] > B_running_max_old_double[,j-1,s],
                                      B_old_double_bottom[,j,s], B_running_max_old_double[,j-1,s])
    Y_old_drawdown[,j,s] <- Y_old_drawdown[,j,s] * pmax(1,B_running_max_old_double[,j,s])
    Y_running_max_old_double[,j,s] <- ifelse(Y_old_drawdown[,j,s] > Y_running_max_old_double[,j-1,s],
                                             Y_old_drawdown[,j,s], Y_running_max_old_double[,j-1,s])

  }



} 


cat("\nFinished.")



#=====================#

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#For the GNU General Public License, see <https://www.gnu.org/licenses/>.
