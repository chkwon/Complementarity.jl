# Complementarity.jl

<!-- [![Build Status](https://travis-ci.org/chkwon/Complementarity.jl.svg?branch=master)](https://travis-ci.org/chkwon/Complementarity.jl) -->


This package provides a modeling and computational interface for solving Mixed Complementarity Problems (MCP): modeling by [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl) and computing by [PATHSolver.jl](https://github.com/chkwon/PATHSolver.jl).


# OS X

At this moment, `PATHSolver.jl` is only available for Mac OS X, and its installation requires some work. Because this package `Complementarity.jl` is dependent on `PATHSolver.jl`, it also is available for Mac OS X only.

# Installation

Currently, this package has not yet been registered.
```julia
Pkg.clone("https://github.com/chkwon/Complementarity.jl.git")
```


# Example

```julia
using Complementarity, JuMP

M = [0  0 -1 -1 ;
     0  0  1 -2 ;
     1 -1  2 -2 ;
     1  2 -2  4 ]

q = [2; 2; -2; -6]

lb = zeros(4)
ub = Inf*ones(4)

m = MCPModel()
@defVar(m, x[1:4])
@defNLExpr(m, F[i=1:4], sum{M[i,j]*x[j], j=1:4} + q[i])
correspond(m, lb, x, ub, F)

PATHSolver.path_options(   
                "convergence_tolerance 100",
                "output no",
                "time_limit 3600"      )

z, f = solveMCP(m)
````
The result should be `z == [2.8, 0.0, 0.8, 1.2]`.

```julia
m = MCPModel()
```
This line prepares a JuMP Model, just same as in [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl).

```julia
@defVar(m, x[1:4])
```
Defining variables is exactly same as in JuMP.jl.

```julia
@defNLExpr(m, F[i=1:4], sum{M[i,j]*x[j], j=1:4} + q[i])
```
This is to define expressions for `F` in MCP. Even when the expression is linear or quadratic, you should use the nonlinear version `@defNLExpr`.

```julia
correspond(m, lb, x, ub, F)
```
This function matches each element of `x` and the corresponding element of `F`.

```julia
PATHSolver.path_options(   
                "convergence_tolerance 100",
                "output no",
                "time_limit 3600"      )
```
This adjusts options of the PATH Solver. See the [list of options](http://www.cs.wisc.edu/~ferris/path/options.pdf).

```
z, f = solveMCP(m)
```
This solves the MCP and receives solution `z` and its function value `f`. I used here `z` and `f`, instead of `x` and `F` not to be confused with the JuMP variable `x` and JuMP expression `F`.
