
# Branch Information

This is a branch of Complementarity. The goal of this branch is to update Complementarity to work with the newest version of JuMP. Currently the only broken test is MCEP.  


# Complementarity.jl


[![Build Status](https://github.com/chkwon/Complementarity.jl/workflows/CI/badge.svg?branch=master)](https://github.com/chkwon/Complementarity.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/chkwon/Complementarity.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/chkwon/Complementarity.jl)


This package provides modeling language for (1) mixed complementarity problems (MCP) and (2) mathematical programs with equilibrium problems (MPEC).

**NOTE** `@complmentarity` for MCP and `@complements` for MPEC.

## Mixed Complementarity Problems (MCP)

**NOTE:** The [PATHSolver.jl](https://github.com/chkwon/PATHSolver.jl) has been completely rewritten between `v0.6.2` and `v1.1.0`. Now PATHSolver.jl provides both an interface to the [PATH](http://pages.cs.wisc.edu/~ferris/path.html) solver and an integration to JuMP, but only limited to *linear* problems at this moment. For *nonlinear* problems, you still need to use [Complementarity.jl](https://github.com/chkwon/Complementarity.jl), which now also uses the new PATHSolver.jl as its solver. Most parts of Complementarity.jl remain the same, except how the solver options are passed. 


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

**NOTE:** For solving MPEC, JuMP.jl `v0.21` has started supporting [complementarity constraints](http://www.juliaopt.org/JuMP.jl/v0.21.1/constraints/#Complementarity-constraints-1). At this moment, [GAMS.jl](https://github.com/GAMS-dev/gams.jl) and [KNITRO](https://github.com/jump-dev/KNITRO.jl) support complementarity constraints.


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
m = Model(Ipopt.Optimizer)
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
