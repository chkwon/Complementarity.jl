**NOTE** `JuMP.jl` `v0.21` has started supporting [complementarity constraints](http://www.juliaopt.org/JuMP.jl/v0.21.1/constraints/#Complementarity-constraints-1). It is unclear which solvers support optimization problems with complementarity constraints (or MPEC) yet. For solving Mixed Complementarity Problems via the PATH solver, watch [PATH.jl](https://github.com/odow/PATH.jl), which provides a Julia native access to the PATH solver and integration with JuMP.jl.


# Complementarity.jl

[![Build Status](https://travis-ci.org/chkwon/Complementarity.jl.svg?branch=master)](https://travis-ci.org/chkwon/Complementarity.jl)
[![Coverage Status](https://coveralls.io/repos/github/chkwon/Complementarity.jl/badge.svg?branch=master)](https://coveralls.io/github/chkwon/Complementarity.jl?branch=master)

This package provides modeling language for (1) mixed complementarity problems (MCP) and (2) mathematical programs with equilibrium problems (MPEC).

**NOTE** `@complmentarity` for MCP and `@complements` for MPEC.

## Mixed Complementarity Problems (MCP)

***[MCP Documentation](MCP.md)***

- This package provides a modeling and computational interface for solving [Mixed Complementarity Problems](https://en.wikipedia.org/wiki/Mixed_complementarity_problem) (MCP): modeling by [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl) and computing by [PATHSolver.jl](https://github.com/chkwon/PATHSolver.jl) and [NLsolve.jl](https://github.com/EconForge/NLsolve.jl). See [the documentation](MCP.md).

```
F(x) ⟂ lb ≤ x ≤ ub
```

A very simple example:
```
(x+2) x = 0,  x ≥ 0,   x+2 ≥ 0
```

```julia
using Complementarity, JuMP
m = MCPModel()
@variable(m, x >= 0)
@mapping(m, F, x+2)
@complementarity(m, F, x)
status = solveMCP(m)
@show result_value(x)
```


## Mathematical Programs with Equilibrium Constraints (MPEC)

***[MPEC Documentation](MPEC.md)***

- For solving [mathematical programs with equilibrium constraints (MPEC)](https://en.wikipedia.org/wiki/Mathematical_programming_with_equilibrium_constraints), this package provides an extension to [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl) by providing a macro that accepts complementarity conditions as constraints.  Then it reformulates the complementarity conditions as a set of equality and inequality constraints so that a nonlinear optimization solver such as [Ipopt.jl](https://github.com/JuliaOpt/Ipopt.jl) can solve the problem. See [the documentation](MPEC.md).

```
min  f(x)
s.t. g(x) ≤ 0
     F(x) ⟂ lb ≤ x ≤ ub
```

A very simple example:
```
min  x^3
s.t. (x+2) x = 0,  x ≥ 0,   x+2 ≥ 0
```

```julia
using JuMP, Ipopt, Complementarity
m = Model(solver=IpoptSolver())
@variable(m, x>=0)
@NLobjective(m, Min, x^3)
@complements(m, 0 <= x+2,   x >= 0)
solve(m)
@show getvalue(x)
```

# Installation

```julia
Pkg.add("Complementarity")
```

This will also install a few other packages.
