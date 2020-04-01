# Using Convex.jl to test if the probability distributions in Homework problem 1
# Sets up the optimization problem for the Cut Inflation for the triangle causal
# model with the distribution:
#
# P(a, b, c) = if (a == b == c) => 1/2
#              else             => 0

using Convex
using GLPK
using Inflation

# the probability matrix is set up as
# P_{A2,B1} = sum_{C1}P_{A2,B1,C1}   --> 'yyx'  (3rd term is integrated)
# P_{A2,C1} = sum_{B1}P_{A2,B1,C1}   --> 'yxy'  (2nd term is integrated)
# P_{B1,C1} = sum_{A2}P_{A2,B1,C1}   --> 'xyy'  (1st term is integrated)
infl_strings = ["yyx", "yxy", "xyy"]

# creating the M-matrix and b-array to be used in the Convex.jl optimization
M = get_M_matrix(infl_strings)
b = get_b_array(infl_strings, pfunc1, cut)

# the variables that get optimized by Convex.jl
v = get_v_array_Variable(infl_strings)

# set up the objective and constraints
objective   = 0
constraints = [v >= 0, M*v == b]

# create the problem and solve
problem = Convex.maximize(objective, constraints)
solve!(problem, GLPK.Optimizer, verbose = false)

# return the results
if problem.status == Convex.MathOptInterface.INFEASIBLE
	println("Solution is infeasible.")
else
	println("Solution is feasible.")
	println("objective = $(problem.optval)")
	println("v = $(v.value)")
end
