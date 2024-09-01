There are 3 files:

(1) Trailing barrier - fractional drawdown via recursive reflection of entire sequence - for github.R

(2) Trailing barrier - fractional drawdown via recursive reflection of each data point - for github.R

(3) Draw trailing barrier - for github.R

----------

(1) was used to produce the sequence of diagrams in my note 11 May 2024 note (subsequently abbreviated to the 2 panels of Figure 7 in my main note.

(2) is more concise code, but equivalent (repeatedly reflecting one data point at a time = repeatedly reflecting the whole sequence.)  

(3) draws the graphs from either. At the top, RunNo <- 33 specifies the simulation correspnding to the graphs in my note (at least on my hardware).  
      - for Figure 6, Z_t_over_M_t is the fractional drawdown F_t, directly-constructed by dividing blue asset path by its own maximum M_t (sorry about the wrong "Z_t" label)  
      - for Figure 7, Y_old_drawdown is the factional drawdwon F_t, constructed by reflecting the red asset path (unrestrited GBM) at 0 and 1. 

The parameter D <- 2 near the top of the main code specifies the number of recursions of the {lower and upper} pair of reflection operators.This needs to be at least: the number of complete dradowns, i.e. number of times the fractional drawdown traverses from 1 to b. So if you set b closer to 1, you will probably need to increase D (because more traversals from 1 to be are likely when 1 and b are close). 

If you set D to a value greater than the number of complete drawdowns, this does no harm: the later applications of reflection operators have no effect, because the process has already been constrained to be between b and 1.  


