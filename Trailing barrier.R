#=====================#
# DISCLAIMER
# This  code is provided "as is" and "with all faults" for educational and research use only.  
#   The author makes no representations or warranties of any kind concerning its correctness or 
#   suitability for any particular purposes.  Some numbers and graphs in published papers were
#   produced in Excel, not this code. If you find errors or ways to improve the code, please 
#   let me know at  R.G.ThomasNOSPAM AT kent.ac.uk
#=====================#

# Setting parameters 

R <- 20 # arbitrary "maximum conceivable number of times drawdown traverse from 1 to b (top to bottom)

S <- 1.0 #stock price at time t
K <- 1.0 #strike price  
b <- 0.5 #barrier - for b = K case, set b just slightly less than K, to avoid divisions by zero problems.
tau <- 10 #time to maturity in years
r <- 0.015 #risk-free annual interest rate 
q <- 0.01 #deferment rate (yield) (needs slightly different from r, to avoid division-by-zero problems in theta)
sigma <- 0.30 #annual volatility of the stock price (standard deviation)

c <- -0.0001 # to investigate limit on trading close to barrier: cannot trade when Y - b < c. 
#It turns out c has only marginal effect
#If don't want limit, set c to tiny neg number is safest (tiny discrepancies otherwise, prob from "Y = b exactly" cases).  
h <- 1 # hedging frequency - hedge every h time steps. Normally I just use 1.

set.seed(1930) #set the seed (if not set, it's taken from the computer clock)
N <- 6300 #N is number of time steps, e.g. 252 x 25 = 6300 steps for every working day over 25 years. 
nSim <- 100 #number of simulations (paths) 


#Check validity of inputs
stopifnot(b <= min(S,K), r!=q, K!=b)

#analytic prices
# z's as in the paper
z1 <- (log(S/K) + (r - q + 0.5*sigma^2)*tau)/(sigma*sqrt(tau))
z2 <- (log(b^2/(K*S)) + (r - q + 0.5*sigma^2)*tau)/(sigma*sqrt(tau))
z3 <- (log(S/b) + (r - q + 0.5*sigma^2)*tau)/(sigma*sqrt(tau))
z4 <- (log(b/S) + (r - q + 0.5*sigma^2)*tau)/(sigma*sqrt(tau))
theta <- 2*(r-q)/sigma^2


BS_call <- S*exp(-q*tau)*pnorm(z1) - K*exp(-r*(tau))*pnorm(z1-sigma*sqrt(tau)) #BS value call
BS_put <- -S*exp(-q*tau)*pnorm(-z1) + K*exp(-r*tau)*pnorm(-z1+sigma*sqrt(tau)) #BS value put

Analytic_barrier_call <- S*exp(-q*tau)*pnorm(z1) - K*exp(-r*(tau))*pnorm(z1-sigma*sqrt(tau))+
1/theta*(S*exp(-q*tau)*(b/S)^(1+theta)*pnorm(z2) - K*exp(-r*(tau))*(K/b)^(theta-1)*pnorm(z2-theta*sigma*sqrt(tau)))

Analytic_barrier_put <- K*exp(-r*tau)*pnorm(-z1+sigma*sqrt(tau)) - S*exp(-q*tau)*pnorm(-z1)-
b*exp(-r*tau)*pnorm(-z3+sigma*sqrt(tau)) + S*exp(-q*tau)*pnorm(-z3)+
1/theta*(b*exp(-r*tau)*pnorm(-z3+sigma*sqrt(tau)) - S*exp(-q*tau)*(b/S)^(1+theta)*(pnorm(z4)-pnorm(z2)) - K*exp(-r*tau)*(K/b)^(theta-1)*pnorm(z2-theta*sigma*sqrt(tau)))


# Compute the Monte Carlo prices 

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

B_double_top <- array(dim=c(nSim,N,R)) # extra s-dimension is to give loops for recursive reflection
B_double_bottom <- array(dim=c(nSim,N,R))
B_running_min_double <- array(dim=c(nSim,N,R))
B_running_max_double <- array(dim=c(nSim,N,R))
X_double <- array(dim=c(nSim,N,R))
Y_running_min_double <- array(dim=c(nSim,N,R))
Y_running_max_double <- array(dim=c(nSim,N,R))

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

X_double[,1,] <- X[,1]
B_running_min_double[,1,] <- 1 /X[,1]
B_running_max_double[,1,] <- b_as_matrix[,1]/X[,1]
Y_running_min_double[,1,] <- S
Y_running_max_double[,1,] <- S

Option_payoff_written_put <- numeric(nSim)

Option_exercise_payment_written_put <- numeric(nSim)


BS_final_cash <-numeric(nSim)
BS_final_stock <-numeric(nSim)
BS_hedging_error <- numeric(nSim)

Thomas_final_cash <-numeric(nSim)
Thomas_final_stock <-numeric(nSim)
Thomas_hedging_error <- numeric(nSim)


#Black-Scholes replication of barrier put    

BS_delta_written_put <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
BS_delta_written_put_most_recent_permitted_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
BS_stock_position_before_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
BS_stock_position_after_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
BS_trade_size <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
BS_cash_account <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
  
BS_delta_written_put[,1] <- -exp(-q*tau)*(pnorm(z1)-1) #positive number
BS_stock_position_after_trade[,1] <- BS_delta_written_put[1] 
BS_cash_account[,1] <-  -S * BS_delta_written_put[1] - BS_put 
BS_delta_written_put_most_recent_permitted_trade[,1] <-BS_delta_written_put[,1] #opening the position always permitted
BS_trade_size[,1] <- 0
  
#Thomas put replication - going to flip the signs at then end to get BOUGHT put
Thomas_delta_written_put <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_delta_written_put_most_recent_permitted_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_stock_position_after_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_trade_size <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_cash_account <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
  
Thomas_z1 <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_z2 <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_z3 <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_z4 <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
  
Thomas_delta_written_put[,1] <- -exp(-q*tau)*(pnorm(z1) - pnorm(z3) + (b_trailing[,1]/S)^(1+theta)*(pnorm(z4)-pnorm(z2)))
Thomas_stock_position_after_trade[,1] <- Thomas_delta_written_put[,1] # same as above
Thomas_cash_account[,1] <- -S * Thomas_delta_written_put[,1] - Analytic_barrier_put #negative number, borrowing to buy the asset. 
Thomas_delta_written_put_most_recent_permitted_trade[,1] <-Thomas_delta_written_put[,1] #opening the position always permitted
Thomas_trade_size[,1] <- 0
  
#Thomas call replication

Thomas_delta_bought_call <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_delta_bought_call_most_recent_permitted_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_call_stock_position_after_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_call_trade_size <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_call_cash_account <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_call_replicating_portfolio <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)

Thomas_delta_bought_call[,1] <- exp(-q*tau)*(pnorm(z1) - (b_trailing[,1]/S)^(1+theta)*pnorm(z2))
#using z1, z2 from top of program...OK at time 1.
Thomas_call_stock_position_after_trade[,1] <- Thomas_delta_bought_call[,1] # same as above
Thomas_call_cash_account[,1] <- -S * Thomas_delta_bought_call[,1] + Analytic_barrier_call 
#Negative number, as for BS, for low barrier. Can be positive for high barrier (when
#the delta of the Thomas call gets low).
Thomas_delta_bought_call_most_recent_permitted_trade[,1] <-Thomas_delta_bought_call[,1] #opening the position always permitted
Thomas_call_trade_size[,1] <- 0



  
#Then do the loop - with hedging every h time steps
  
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
    

   
      if(j %% h > 0) {
        BS_cash_account[,j] <- BS_cash_account[,j-1] * exp(r*dt) + BS_stock_position_after_trade[,j-1] *Y[,j] * q *dt 
        BS_stock_position_after_trade[,j] <- BS_stock_position_after_trade[,j-1]
        BS_delta_written_put_most_recent_permitted_trade[,j] <- BS_delta_written_put_most_recent_permitted_trade[,j-1]
        
      }else{
        
        BS_delta_written_put[,j] <- -exp(-q*tau*(1-j/N))*(pnorm((log(Y[,j]/K) +(r- q +0.5*sigma^2)*(tau*(1-j/N)))/(sigma*sqrt(tau*(1-j/N))))-1)
        
        BS_trade_size[,j] <- ifelse(Y[,j] - b_trailing[,j] < c, 0, BS_delta_written_put[,j] - BS_delta_written_put_most_recent_permitted_trade[,j-1])
        BS_stock_position_after_trade[,j] <- BS_stock_position_after_trade[,j-1] + BS_trade_size[,j] 
        BS_cash_account[,j] <- BS_cash_account[,j-1] * exp(r*dt) + BS_stock_position_after_trade[,j-1] *Y[,j] * q *dt - BS_trade_size[,j] * Y[,j] 
        BS_delta_written_put_most_recent_permitted_trade[,j] <- ifelse(Y[,j] - b < c, BS_delta_written_put_most_recent_permitted_trade[,j-1],
              BS_delta_written_put[,j]) 
      }    


    #Then for Thomas put hedging
    
    Thomas_z1[,j] <-(log(Y[,j]/K) +(r - q +0.5*sigma^2)*(tau*(1-j/N)))/(sigma*sqrt(tau*(1-j/N)))
    Thomas_z2[,j] <-(log(b_trailing[,j]^2/(K*Y[,j])) +(r - q +0.5*sigma^2)*(tau*(1-j/N)))/(sigma*sqrt(tau*(1-j/N)))
    Thomas_z3[,j] <-(log(Y[,j]/b_trailing[,j]) +(r - q + 0.5*sigma^2)*(tau*(1-j/N)))/(sigma*sqrt(tau*(1-j/N)))
    Thomas_z4[,j] <-(log(b_trailing[,j]/Y[,j]) +(r - q + 0.5*sigma^2)*(tau*(1-j/N)))/(sigma*sqrt(tau*(1-j/N)))
 
    if(j %% h > 0) {
      Thomas_cash_account[,j] <- Thomas_cash_account[,j-1] * exp(r*dt) + Thomas_stock_position_after_trade[,j-1] *Y[,j] * q *dt
      Thomas_stock_position_after_trade[,j] <- Thomas_stock_position_after_trade[,j-1]
      Thomas_delta_written_put_most_recent_permitted_trade[,j] <- Thomas_delta_written_put_most_recent_permitted_trade[,j-1]
    
       }else{
    
           if (j!=N) {
           Thomas_delta_written_put[,j] <- -exp(-q*tau*(1-j/N))*(pnorm(Thomas_z1[,j]) - pnorm(Thomas_z3[,j]) +
             (b_trailing[,j]/Y[,j])^(1+theta) * (pnorm(Thomas_z4[,j]) - pnorm(Thomas_z2[,j])))
             } else {
           Thomas_delta_written_put[,j] <- ifelse(K > Y[,j], 1, 0)
            }
    
       Thomas_trade_size[,j] <- ifelse(Y[,j] - b_trailing[,j] < c, 0, Thomas_delta_written_put[,j] - Thomas_delta_written_put_most_recent_permitted_trade[,j-1])
       Thomas_stock_position_after_trade[,j] <- Thomas_stock_position_after_trade[,j-1] + Thomas_trade_size[,j]
       Thomas_cash_account[,j] <- Thomas_cash_account[,j-1] * exp(r*dt) + Thomas_stock_position_after_trade[,j-1] *Y[,j] * q *dt -
                               Thomas_trade_size[,j] * Y[,j]
    
       Thomas_delta_written_put_most_recent_permitted_trade[,j] <- ifelse(Y[,j] - b_trailing[,j] < c, 
             Thomas_delta_written_put_most_recent_permitted_trade[,j-1], Thomas_delta_written_put[,j])
    
       }

    
     #Then for Thomas call hedging 
    
    # Already calculated the required column of z_1[,j] and column of z_2[,j] above

    if(j %% h > 0) {
      Thomas_call_cash_account[,j] <- Thomas_call_cash_account[,j-1] * exp(r*dt) + Thomas_call_stock_position_after_trade[,j-1] *Y[,j] * q *dt
      Thomas_call_stock_position_after_trade[,j] <- Thomas_call_stock_position_after_trade[,j-1]
      Thomas_delta_bought_call_most_recent_permitted_trade[,j] <- Thomas_delta_bought_call_most_recent_permitted_trade[,j-1]    
     
       }else{
      
      
          if (j!=N) {
          Thomas_delta_bought_call[,j] <- exp(-q*tau*(1-j/N))*(pnorm(Thomas_z1[,j]) -
                                                             (b_trailing[,j]/Y[,j])^(1+theta) * pnorm(Thomas_z2[,j]))
           } else {
          Thomas_delta_bought_call[,j] <- ifelse(K < Y[,j], 1, 0)
          }
    # Setting the final delta manually avoids the deep bug of z2, z3, z4 evaluating as 0/0 if Y_T[i,j]] = b_trailing[,j] 
    #   and time remaining = 0. 
    
       
         Thomas_call_trade_size[,j] <- ifelse(Y[,j] - b_trailing[,j] < c, 0, Thomas_delta_bought_call[,j] - Thomas_delta_bought_call_most_recent_permitted_trade[,j-1])
         Thomas_call_stock_position_after_trade[,j] <- Thomas_call_stock_position_after_trade[,j-1] + Thomas_call_trade_size[,j]
         Thomas_call_cash_account[,j] <- Thomas_call_cash_account[,j-1] * exp(r*dt) + Thomas_call_stock_position_after_trade[,j-1] *Y[,j] * q *dt -
                                          Thomas_call_trade_size[,j] * Y[,j]
         
         Thomas_delta_bought_call_most_recent_permitted_trade[,j] <- ifelse(Y[,j] - b_trailing[,j] < c, 
                 Thomas_delta_bought_call_most_recent_permitted_trade[,j-1], Thomas_delta_bought_call[,j])        
        
            
        }
    
    
  

#=================#
   
}  # end of N steps

# Now construct the drawdown of the "GBM reflected at trailing barrier" by recursive reflection 
# of the original GBM at 1 and b. 
# NO problem if (exceptionally) GBM is already confined between b an 1 before the end of the 
#minimum of two (upper, lower) applications of reflection below -in that case, the further reflections
# just have no effect.

for (j in 2:N){ # Reflect original GBM, X, at the top first:

B_double_top[,j,1] <- 1 /X[,j]
B_running_min_double[,j,1] <- ifelse(B_double_top[,j,1] < B_running_min_double[,j-1,1],
                                         B_double_top[,j,1], B_running_min_double[,j-1,1])
X_double[,j,1] <- X[,j] * pmin(1,B_running_min_double[,j,1])
Y_running_min_double[,j,1] <- ifelse(X_double[,j,1] < Y_running_min_double[,j-1,1],
                                         X_double[,j,1], Y_running_min_double[,j-1,1])

}

for (j in 2:N){ # Then reflect result at the bottom:

  B_double_bottom[,j,1] <- b_as_matrix[,j]/X_double[,j,1] # from above
  B_running_max_double[,j,1] <- ifelse(B_double_bottom[,j,1] > B_running_max_double[,j-1,1],
                                           B_double_bottom[,j,1], B_running_max_double[,j-1,1])
  X_double[,j,1] <- X_double[,j,1] * pmax(1,B_running_max_double[,j,1])
  Y_running_max_double[,j,1] <- ifelse(X_double[,j,1] > Y_running_max_double[,j-1,1],
                                           X_double[,j,1], Y_running_max_double[,j-1,1])

}

# Now reflect the result at {top, bottom} for up another R-1 times
# (R is the maximum allowed number of times the drawdown traverses 1 to b -  in practice break below when
# all paths in the simulation are OK.)


for (s in 2:R){

   for (j in 2:N){   # Reflect original GBM, X, at the top first:

    B_double_top[,j,s] <- 1 /X_double[,j,s-1]  # has to be s-1 here, not found for s yet
    B_running_min_double[,j,s] <- ifelse(B_double_top[,j,s] < B_running_min_double[,j-1,s],
                                             B_double_top[,j,s], B_running_min_double[,j-1,s])
    X_double[,j,s] <- X_double[,j,s-1] * pmin(1,B_running_min_double[,j,s])
    Y_running_min_double[,j,s] <- ifelse(X_double[,j,s] < Y_running_min_double[,j-1,s],
                                             X_double[,j,s], Y_running_min_double[,j-1,s])

  }

  for (j in 2:N){ # Then reflect result at the bottom (and so on, recursively)
    
    B_double_bottom[,j,s] <- b_as_matrix[,j]/X_double[,j,s] # s as final subscript here, using output from directly above
    B_running_max_double[,j,s] <- ifelse(B_double_bottom[,j,s] > B_running_max_double[,j-1,s],
                                             B_double_bottom[,j,s], B_running_max_double[,j-1,s])
    X_double[,j,s] <- X_double[,j,s] * pmax(1,B_running_max_double[,j,s])
    Y_running_max_double[,j,s] <- ifelse(X_double[,j,s] > Y_running_max_double[,j-1,s],
                                             X_double[,j,s], Y_running_max_double[,j-1,s])
    
  }

Job_done <- ifelse(max(X_double[,,s]) - min(X_double[,,s]) > 1-b, 0,1)
  
if(max(Job_done) > 0){break} # job is done when all paths of X_double are confined between 1 and b

} #end of s-loop


Y_T <- Y[,N]
  
 
Option_payoff_written_put <- - pmax(K - Y_T,0) # minus sign as _written_ put. pmax is "parallel maximum" of 2 vectors.
Option_payoff_bought_call <- pmax(Y_T - K,0)
  
#First for BS replication

BS_final_cash <- BS_cash_account[,j]
BS_final_stock <-  BS_stock_position_after_trade[,j]


BS_hedging_error <- ((BS_final_cash + BS_final_stock * Y_T) - Option_payoff_written_put)
# ("proceeds of replication scheme" - "option payoff")

#Then for Thomas written put replication
  # NB NB Going to FLIP SIGNS of final stock and final cash here, because
  # we want to get payoffs of BOUGHT put on the Thomas side

Thomas_final_cash <- (-1) * Thomas_cash_account[,j]

Thomas_final_stock <- (-1) * Thomas_stock_position_after_trade[,j]

Thomas_hedging_error <- ((Thomas_final_cash + Thomas_final_stock * Y_T) + Option_payoff_written_put)
  # flipped sign on option payoff here, because we want payoff of a BOUGHT put here, having
  # converted to the payoffs of a BOUGHT put by x (-1) on stock and cash just above.

#Then for Thomas bought call replication

Thomas_call_final_cash <- Thomas_call_cash_account[,j]
Thomas_call_final_stock <-  Thomas_call_stock_position_after_trade[,j]


Thomas_call_hedging_error <- ((Thomas_call_final_cash + Thomas_call_final_stock * Y_T) - Option_payoff_bought_call)
  


payoff_expiry_call <-pmax(Y_T-K,0) 
expected_payoff_call <- mean(payoff_expiry_call)
Monte_Carlo_call_price <- exp(-r*(tau))*expected_payoff_call

payoff_expiry_put <-pmax(K-Y_T,0) 
expected_payoff_put <- mean(payoff_expiry_put)
Monte_Carlo_put_price <- exp(-r*(tau))*expected_payoff_put

#First meansn for BS 
Mean_BS_final_cash <-mean(BS_final_cash)
Max_BS_final_cash <- max(BS_final_cash)
Min_BS_final_cash <- min(BS_final_cash)
Mean_BS_hedging_error <- mean(abs(BS_hedging_error))

#Then means for Thomas put
Thomas_mean_final_cash <-mean(Thomas_final_cash)
Thomas_max_final_cash <- max(Thomas_final_cash)
Thomas_min_final_cash <- min(Thomas_final_cash)
Thomas_mean_hedging_error <- mean(abs(Thomas_hedging_error))

#Then means for Thomas call
Thomas_call_mean_final_cash <-mean(Thomas_call_final_cash)
Thomas_call_max_final_cash <- max(Thomas_call_final_cash)
Thomas_call_min_final_cash <- min(Thomas_call_final_cash)
Thomas_call_mean_hedging_error <- mean(abs(Thomas_call_hedging_error))

Ratio_call <- Monte_Carlo_call_price/Analytic_barrier_call
Ratio_put <- Monte_Carlo_put_price/Analytic_barrier_put

mybins <-seq(0.0, 0.1,0.001) #This can be used to make both histograms look same, with narrow bins.
mybins2 <-seq(-0.9, 0.1, 0.01) 
mybins3 <-seq(-0.6, 0, 0.01)
hist(BS_hedging_error)
hist(Thomas_hedging_error)
hist(ND_final_cash) #, breaks=mybins)
hist(ND_lowest_cash) # , breaks=mybins2)


#Calculate the theoretical gains from the trading - exactly offsetting the rolled-forward difference in initial valuations
Rolled_up_difference_in_call_valuations <- (Analytic_barrier_call - (S * exp(-q*tau) - K*exp(-r*tau) + Analytic_barrier_put)) * exp(r * tau)

cat("Analytic B-Scholes Call :",BS_call, "   Analytic B-Scholes Put :",BS_put)

cat("\nAnalytic Barrier Call   :",Analytic_barrier_call, "  Analytic Barrier Put   :",Analytic_barrier_put)

cat("\nMonte-Carlo Barrier Call:",Monte_Carlo_call_price, " Monte-Carlo Barrier Put:",Monte_Carlo_put_price)

cat("\nMean B-Scholes put abs hedging error:",round(Mean_BS_hedging_error,digits=7),  "Thomas ditto: ",round(Thomas_mean_hedging_error,digits=7))

cat("\nMean Thomas call abs hedging error:",round(Thomas_call_mean_hedging_error,digits=7))


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
