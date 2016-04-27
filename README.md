# Complementarity.jl

<!-- [![Build Status](https://travis-ci.org/chkwon/Complementarity.jl.svg?branch=master)](https://travis-ci.org/chkwon/Complementarity.jl) -->


This package provides a modeling and computational interface for solving Mixed Complementarity Problems (MCP): modeling by [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl) and computing by [PATHSolver.jl](https://github.com/chkwon/PATHSolver.jl).


# OS X

At this moment, `PATHSolver.jl` is only available for Mac OS X, and its installation requires some work. Because this package `Complementarity.jl` is dependent on `PATHSolver.jl`, it also is available for Mac OS X only.

# installation
1. Install `PATHSolver.jl`. Read the [instruction of PATHSolver.jl](https://github.com/chkwon/PATHSolver.jl).
2. Install `Complementarity.jl`:
```julia
Pkg.clone("https://github.com/chkwon/Complementarity.jl.git")
```

# Example

```julia
using Complementarity, JuMP

m = MCPModel()

M = [0  0 -1 -1 ;
     0  0  1 -2 ;
     1 -1  2 -2 ;
     1  2 -2  4 ]

q = [2; 2; -2; -6]

lb = zeros(4)
ub = Inf*ones(4)

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
