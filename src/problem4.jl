# Using Convex.jl to test if the probability distributions in Homework problem 1
# Sets up the optimization problem for the Cut Inflation for the triangle causal
# model with the distribution:
#
# P(a, b, c) = if (a + b + c == 1) => 1/3
#              else             => 0

using Convex
using GLPK
using Inflation

# the probability matrix is set up as
# P_{A1,B1,C1} = sum_{A2,B2,C2}P     --> 'yxyxyx'
# P_{A1,B2,C2} = sum_{A2,B1,C1}P     --> 'yxxyxy'
# P_{A2,B1,C2} = sum_{A1,B2,C1}P     --> 'xyyxxy'
# P_{A2,B2,C1} = sum_{A1,B1,C2}P     --> 'xyxyyx'
# P_{A2,B2,C2} = sum_{A1,B1,C1}P     --> 'xyxyxy'
infl_strings = ["yxyxyx", "yxxyxy", "xyyxxy", "xyxyyx", "xyxyxy"]

# creating the M-matrix and b-array to be used in the Convex.jl optimization
M = get_M_matrix(infl_strings)
b = get_b_array(infl_strings, pfunc2, spiral)

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
	println("Solution is infeasible.")
else
	println(problem.optval)
	println(v.value)
end
