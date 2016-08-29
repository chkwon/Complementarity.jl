
# Mathematical Programming with Equilibrium Constraints (MPEC)


The main objective of this package is to provide the `complements` keyword (macro in this Julia package) that is available in other modeling languages such as GAMS, Pyomo, and AMPL.


# `@complements` macro

```julia
using JuMP, Ipopt, Complementarity
m = Model(solver=IpoptSolver())
```

The complementarity must be written one of the following forms:
```julia
@complements(m, F,   lb <= x <= ub)
@complements(m, lb <= F,        x >= lb)
@complements(m, ub >= F,        x >= lb)
@complements(m, lb <= F,        x <= ub)
@complements(m, ub >= F,        x <= ub)
```
where `F` is a JuMP expression and `x` is a JuMP variable. (**NOTE:** Please be careful with the location of the inequalities; bounds on `F` needs to be on the left-hand side, and bounds on `x` needs to be on the right-hand-side.)


In most cases, `lb` is zero for both `F` and `x`, so likely:
```julia
@complements(m, 0 <= F,        x >= 0)
```

Vector or block inputs to `@complements` are not supported yet. That is, you have to do something like
```julia
for i in I
    @complements(m, 0 <= F[i],        x[i] >= 0)
end
```
instead of
```julia
# This code does NOT work. It's in the wish list.
@complements(m, cc[i=1:n], 0 <= F[i],   x[i] >= 0)
```

Then, the MPEC.jl package will transform the complementarity conditions provided to a set of nonlinear equality and inequality constraints.

See the examples below.

# Examples

### `bard1.jl` --- translated from `bard1.mod` in [MacMPEC](https://wiki.mcs.anl.gov/leyffer/index.php/MacMPEC)
```julia
using JuMP, Ipopt, Complementarity

m = Model(solver=IpoptSolver())
@variable(m, x>=0)
@variable(m, y>=0)
@variable(m, l[1:3])

@NLobjective(m, Min, (x - 5)^2 + (2*y + 1)^2)

@NLconstraint(m, 2*(y-1) - 1.5*x + l[1] - l[2]*0.5 + l[3] == 0)

@complements(m, 0 <= 3*x - y - 3,        l[1] >= 0)
@complements(m, 0 <= - x + 0.5*y + 4,    l[2] >= 0)
@complements(m, 0 <= - x - y + 7,        l[3] >= 0)

solve(m)

@show getobjectivevalue(m)
@assert isapprox(getobjectivevalue(m), 17.0000, atol=1e-4)
```

### `dempe.jl` --- translated from `dempe.mod` in [MacMPEC](https://wiki.mcs.anl.gov/leyffer/index.php/MacMPEC)
```julia
using JuMP, Ipopt, Complementarity

@mpec_tolerance 100.0

m = Model(solver=IpoptSolver())
@variable(m, x)
@variable(m, z)
@variable(m, w>=0)

@NLobjective(m, Min, (x - 3.5)^2 + (z + 4)^2)

@NLconstraint(m, z - 3 + 2*z*w == 0)

@complements(m, 0 >= z^2 - x,  w >= 0)

solve(m)

@show getobjectivevalue(m)
@assert isapprox(getobjectivevalue(m), 28.25, atol=1e-4)
```


### `gnash1m.jl` --- translated from `gnash1m.mod` in [MacMPEC](https://wiki.mcs.anl.gov/leyffer/index.php/MacMPEC)
```julia
using JuMP, Ipopt, Complementarity

c = [10, 8, 6, 4, 2]
K = [5, 5, 5, 5, 5]
b = [1.2, 1.1, 1.0, 0.9, 0.8]
L = 20
g = 1.7

gg = 5000^(1/g)

m = Model(solver=IpoptSolver())

@variable(m, 0 <= x <= L)
@variable(m, y[1:4])
@variable(m, l[1:4])
@variable(m, Q >= 0)
@constraint(m, Q == x+y[1]+y[2]+y[3]+y[4])
@NLobjective(m, Min, c[1]*x + b[1]/(b[1]+1)*K[1]^(-1/b[1])*x^((1+b[1])/b[1])
 		             - x*( gg*Q^(-1/g) ) )

@NLconstraint(m, cnstr[i=1:4], 0 == ( c[i+1] + K[i+1]^(-1/b[i+1])*y[i] ) - ( gg*Q^(-1/g) )
                               - y[i]*( -1/g*gg*Q^(-1-1/g) ) - l[i] )

for i in 1:4
    @complements(m, l[i], 0 <= y[i] <= L, smooth)
end

solve(m)

@show getobjectivevalue(m)
@assert isapprox(getobjectivevalue(m), -6.11671, atol=1e-4)
```

# Transformations

This package aims to implement transformation techniques mentioned in
- R. Fourer, ["Modeling and Solving Nontraditional Optimization Problems
Session 2b: Complementarity Conditions"](http://ampl.com/MEETINGS/TALKS/2011_01_Chiang_Mai_Workshop_2b.pdf)
- W. E. Hart and J. D. Siirola, ["Modeling Mathematical Programs with Equilibrium Constraints in Pyomo"](http://prod.sandia.gov/techlib/access-control.cgi/2015/155584.pdf)
- M. C. Ferris and S. P. Dirkse, ["Mathematical Programs with Equilibrium Constraints: Automatic Reformulation and
Solution via Constrained Optimization"](http://pages.cs.wisc.edu/%7Eferris/papers/NA-02-11.ps)

At this moment, this package only provides the following simple nonlinear transformation. Other more advanced nonlinear transformation explained in Ferris and Dirkse, and the simple disjunctive transformation described in Hart and Siirola.

## Simple Nonlinear Transformation

For example,
```julia
@complements(m, 0 <= F,        x >= 0, simple)
```
will add the following macros
```julia
@NLconstraint(m, F*x <= mpec_tol)
@NLconstraint(m, x >= 0)
@NLconstraint(m, F >= 0)
```


## Smooth Nonlinear Transformation

This smooth nonlinear transformation is also called Fischer-Burmeister.

For example,
```julia
@complements(m, 0 <= F,        x >= 0, smooth)
```
will add the following macros
```julia
@NLconstraint(m, sqrt( F^2+x^2+mpec_tol) - F - x == 0)
@NLconstraint(m, x >= 0)
@NLconstraint(m, F >= 0)
```

When it is not specified, the default transformation method is `smooth`.


# Some Useful Links

- http://ampl.com/MEETINGS/TALKS/2011_01_Chiang_Mai_Workshop_2b.pdf

- http://www.gamsworld.org/mpec/index.htm
- https://wiki.mcs.anl.gov/leyffer/index.php/MacMPEC
- http://www.neos-guide.org/content/mathematical-programs-equilibrium-constraints
- https://optimization.mccormick.northwestern.edu/index.php/Mathematical_programming_with_equilibrium_constraints#Complementarity_constrained_optimization
- http://www.gamsworld.org/mpec/index.htm
- https://wiki.mcs.anl.gov/leyffer/index.php/MacMPEC
- http://ampl.com/BOOK/CHAPTERS/22-complement.pdf
- http://prod.sandia.gov/techlib/access-control.cgi/2015/155584.pdf
- http://pages.cs.wisc.edu/~ferris/papers/NA-02-11.ps
