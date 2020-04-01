# A simple example to test out Convex.jl

using Convex
using GLPK
using Printf

A = [1.0 1.0 1.0]
c = [1.0, 2.0, -1.0]
b = [1.0]

# the variables that get optimized by Convex.jl
y = Variable(1)
s = Variable(3)

# set up the objective (thing to optimize)
objective   = dot(b, y)

# set up the constrains (things that the Variable instances must obey)
constraints = [s >= 0, s == c - transpose(A)*y]

# create the problem and solve it.
problem = Convex.maximize(objective, constraints)
solve!(problem, GLPK.Optimizer)

# get the results
println(problem.optval)
println(y.value)
println(s.value)
