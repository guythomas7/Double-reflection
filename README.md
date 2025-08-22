There are 2 files:

(1) Trailing barrier - for github.R

(2) Draw trailing barrier - for github.R

----------


The parameter D <- 20 near the top of the main code specifies the number of recursions of the {lower and upper} pair of reflection operators.This needs to be at least: the number of complete drawdowns, i.e. number of times the fractional drawdown traverses from 1 to b. So if you set b closer to 1, you will probably need to increase D (because more traversals from 1 to be are likely when 1 and b are close). 

If you set D to a value greater than the number of complete drawdowns, this does no harm: the later applications of reflection operators have no effect, because the process has already been constrained to be between b and 1.  


