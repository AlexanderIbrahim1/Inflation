# Using Convex.jl to test if the probability distributions in Homework problem 2
# Sets up the optimization problem for the Cut Inflation for the triangle causal
# model with the distribution:
#
# P(a, b, c; t) = if (a == b == c) => t/2
#                 else             => (1 - t)/6
#
# the optimal value of t is found through repetition

using Convex
using GLPK
using Inflation
using Printf

struct InflationLoopData
	all_feasible::Bool
	all_infeasible::Bool
	tmin::Float64
	tmax::Float64
	Ndiv::Int64
end

function check_inflation_t(infl_strings::Vector{String}, M::Array{Int64, 2}, idata::InflationLoopData)::InflationLoopData
	# checks if an inflation is feasible for pfunct, for certain values of t
	tmin = idata.tmin
	tmax = idata.tmax
	Ndiv = idata.Ndiv
	dt   = (tmax - tmin)/(Ndiv - 1)
	for (i, t) in enumerate(range(tmin, tmax, length = Ndiv))
		# avoid repeatedly testing the 1st element
		if i == 1 && t != 0.0
			continue
		end
		
		# get pfunc(a, b, c; t) and get the associated b-array
		pfunct = get_pfunct(t)
		b      = get_b_array(infl_strings, pfunct, cut)
		
		# the variables that get optimized by Convex.jl
		v = get_v_array_Variable(infl_strings)
		
		# set up the objective and constraints
		objective   = 0
		constraints = [v >= 0, M*v == b]
		
		# create the problem and solve
		problem = Convex.maximize(objective, constraints)
		solve!(problem, GLPK.Optimizer, verbose=false)
		
		# return the results
		if problem.status == Convex.MathOptInterface.INFEASIBLE
			@printf("Checked t = %0.6f : solution is infeasible.\n", t)
			if (i == 1) # all infeasible
				return InflationLoopData(false, true, tmin, tmax, Ndiv)
			else
				new_tmin = t - (tmax - tmin)/(Ndiv - 1)
				new_tmax = t
				return InflationLoopData(false, false, new_tmin, new_tmax, Ndiv)
			end
		else
			@printf("Checked t = %0.6f : solution is feasible.\n", t)
		end
	end
	
	return InflationLoopData(true, false, tmin, tmax, Ndiv)
end

# the probability matrix is set up as
# P_{A2,B1} = sum_{C1}P_{A2,B1,C1}   --> 'yyx'  (3rd term is integrated)
# P_{A2,C1} = sum_{B1}P_{A2,B1,C1}   --> 'yxy'  (2nd term is integrated)
# P_{B1,C1} = sum_{A2}P_{A2,B1,C1}   --> 'xyy'  (1st term is integrated)
infl_strings = ["yyx", "yxy", "xyy"]

# creating the M-matrix to be used in the Convex.jl optimization
M = get_M_matrix(infl_strings)

tmin = 0.0  # smallest allowed value of t
tmax = 1.0  # greatest allowed value of t
Ndiv = 11   # number of elements to divide the range into
idata = InflationLoopData(false, false, tmin, tmax, Ndiv)

# check to 5 decimal places
for i in 1:6
	global idata
	idata = check_inflation_t(infl_strings, M, idata)
	if idata.all_feasible
		@printf("All values of t between %0.6f and %0.6f are feasible.\n", tmin, idata.tmax)
		break
	elseif idata.all_infeasible
		@printf("Unable to find further feasible values.\n")
		break
	else
		idata = InflationLoopData(false, false, idata.tmin, idata.tmax, idata.Ndiv)
	end
end

@printf("All values of t between %0.6f and %0.6f are feasible.\n", tmin, idata.tmin)
