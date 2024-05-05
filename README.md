Run the main code.

The parameter R <- 2 near the top specifies the number of recursions of the {lower and upper} pair of reflection operators.This needs to be at least: the number of times the fractional drawdown traverses from 1 to b. So if you set b closer to 1, you will probably need to increase R (because more traversals from 1 to be are likely when 1 and b are close). 

If you set R to a value greater than the number of 1 -> b traversals. this does no harm: the later reflections have no effect, because the process has already been constrained to be between b and 1.  

Then graph the output using the graph code. At the top, RunNo <- 33 specifies the simulation correspnding to the graphs in my note (at least on my hardware).  

