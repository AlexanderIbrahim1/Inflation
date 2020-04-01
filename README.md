# Inflation

## Homework 4 for PHYS 776

A collection of functions that use `Convex.jl` to check the compatibility of probability distributions with different inflations of the triangle causal model.

A probability distribution of the triangle causal model can be studied using linear programming by creating an "inflation" of the model. By integrating different variables in the model, the probability distributions of the inflated model can be written as linear combinations of marginal distributions of the original model. A matrix *M* of *0*'s and *1*'s keeps track of which terms in the linear system are present, and a vector *b* gives the probabilities in terms of the original distribution from the triangle causal model. The goal is to find a vector *v* that satisfies the matrix equation with *M* to give *b*. The elements of *v* represent probabilities, and thus must be non-negative. The problem of finding *v* can be expressed as the following optimization problem (taken from [1]).

![Alt text](./images/optim1.svg)

If the probability distribution applied to an inflation of a causal model is infeasible, that means the original causal model is also infeasible.

## Instructions

The linear equations are set up using strings of *x*'s and *y*'s to represent which terms in a set of equations are integrated and which are evaluated, respectively. For example, with the Cut inflation, we have the equation

P<sub>*A<sub>2</sub>,B<sub>1</sub>*</sub>(*a<sub>2</sub>, b<sub>1</sub>*) = &Sigma;<sub>*c*</sub> P<sub>*ABC*</sub>(*a,b,c*)  .

On the right-hand side, the parameter *c* is integrated over, while the parameters *a* and *b* are evaluated. Thus this set of equations is represented using "abc" => "yyx". The same is done for the other linear sets of equations.

The *M* matrix is constructed from the list of these strings alone. Most of the matrix is *0*'s, and the locations of the *1*'s can be found from binary combinations of the associated strings. For example, "00x" corresponds to the first row in the *M* matrix, and replacing the "x" with "0" and "1" gives "000" and "001", indicating that the *1<sup>st</sup>* and *2<sup>nd</sup>* terms in the first row are *1*'s. Similarly, the *b* vector can be constructed from the same strings, along with the probability distribution function.

## Example

```Julia
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
```

## Problem 1

This problem can be run from `src/problem1.jl`, and checks the feasibility of equation (1) in [1] applied to the Cut inflation for the triangle model. As expected, the distribution and inflation are incompatible.

## Problem 2

This problem can be run from `src/problem2.jl`, and checks the feasibility of equation (16) in [1] applied to the Cut inflation for the triangle model. Different values of *t* between *0* and *1* are used, and by gradually improving the resolution of *t*, it was found that, up to six decimals, distributions using values of *t* between *0.000000* and *0.625000* (or 5/8) are compatible with the Cut inflation.

## Problem 3

This problem can be run from `src/problem3.jl`, and checks the feasibility of equation (17) in [1] applied to the Cut inflation for the triangle model. As expected, the distribution and inflation are compatible. However, the distribution is not compatible with the original triangle causal model [2].

As mentioned in [1], if the linear program is infeasible, it proves that the original system is infeasible. This does not mean that, if the linear program is feasible, that the original system is feasible. As we will see in problem 4, the distribution in equation (17) applied to the Spiral inflation of the triangle model is infeasible.

## Problem 4

This problem can be run from `src/problem4.jl`, and checks the feasibility of equation (17) in [1] applied to the Spiral inflation for the triangle model. As expected, the distribution and inflation are incompatible.

## References
[1] https://github.com/PerimeterInstitute/Computational-Physics-Course-Winter-2020/blob/master/class-2020-03-11/Step2_CausalStructures.pdf
[2] https://github.com/PerimeterInstitute/Computational-Physics-Course-Winter-2020/blob/master/class-2020-03-11/CausalInferencePaper.pdf
