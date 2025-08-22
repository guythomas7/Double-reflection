There are 2 files:

(1) Trailing barrier - for github.R

(2) Draw trailing barrier - for github.R

----------


The parameter D <- 20 near the top of the main code specifies the number of recursions of the {lower and upper} pair of reflection operators.This needs to be at least: the number of complete drawdowns, i.e. number of times the fractional drawdown traverses from 1 to b. The way to think of this is that if the particle is at the latter part of its path, it needs to have received and "remembered" all the "pushes" at the upper and lower barriers earlier in its path.  

If you set D to a value greater than the number of complete drawdowns (e.g 20 as above), this does no harm: the later applications of reflection operators have no effect, because the process has already been constrained to be between b and 1.  


