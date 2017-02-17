# Mixed Complementarity Problems (MCP)

Note that MCP is more general than [Linear Complementarity Problems](https://en.wikipedia.org/wiki/Linear_complementarity_problem) (LCP) and [Nonlinear Complementarity Problems](https://en.wikipedia.org/wiki/Nonlinear_complementarity_problem) (NCP).

The form of MCP is as follows:
```
F(x) ⟂ lb ≤ x ≤ ub
```
which means
- `x = lb`, then `F(x) ≥ 0`
- `lb < x < ub`, then `F(x) = 0`
- `x = ub`, then `F(x) ≤ 0`

When there is no upper bound `ub`, and the lower bound `lb=0`, then it is a regular Nonlinear Complementarity Problem (NCP) of the form:
```
0 ≤ F(x) ⟂ x ≥ 0
```
which means
```
F(x)' x = 0, F(x) ≥ 0, x ≥ 0
```
When `F(x)` is a linear operator such as `F(x) = M x + q` with matrix `M` and vector `q`, then it is a Linear Complementarity Problem (LCP). All these problems are solved by the [PATH Solver](http://pages.cs.wisc.edu/%7Eferris/path.html) which is wrapped by the [PATHSolver.jl](https://github.com/chkwon/PATHSolver.jl) package.

This package `Complementarity.jl` extends the modeling language from [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl) to model complementarity problems.


# Solution via PATHSolver.jl

## License

Without a license, the PATH Solver can solve problem instances up to with up to 300 variables and 2000 non-zeros. For information regarding license, visit the [PATHSolver.jl](https://github.com/chkwon/PATHSolver.jl) page and the [license page](http://pages.cs.wisc.edu/~ferris/path/LICENSE) of the PATH Solver.


## Example 1

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

items = 1:4

# @variable(m, lb[i] <= x[i in items] <= ub[i])
@variable(m, x[i in items] >= 0)
@operator(m, F[i in items], sum{M[i,j]*x[j], j in items} + q[i])
@complementarity(m, F, x)

PATHSolver.options(convergence_tolerance=1e-8, output=:yes, time_limit=3600)


status = solveMCP(m, solver=:PATHSolver)

z = getvalue(x)
````
The result should be `[2.8, 0.0, 0.8, 1.2]`.

```julia
m = MCPModel()
```
This line prepares a JuMP Model, just same as in [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl).

```julia
@variable(m, x[i in items] >= 0)
```
Defining variables is exactly same as in JuMP.jl. Lower and upper bounds on the variables in the MCP model should be provided here.

```julia
@operator(m, F[i in items], sum{M[i,j]*x[j], j in items} + q[i])
```
This is to define expressions for `F` in MCP. This is merely an alias of `JuMP.@NLexpression`.

```julia
@complementarity(m, F, x)
```
This macro matches each element of `F` and the complementing element of `x`.

```julia
PATHSolver.options(convergence_tolerance=1e-8, output=:yes, time_limit=3600)
```
This adjusts options of the PATH Solver. See the [list of options](http://www.cs.wisc.edu/~ferris/path/options.pdf).

```julia
solveMCP(m)
```
This solves the MCP and stores the solution inside `m`, which can be accessed by `getvalue(x)` as in JuMP.


## Example 2

This is a translation of [`transmcp.gms`](http://www.gams.com/modlib/libhtml/transmcp.htm) originally written in GAMS.

```julia
using Complementarity, JuMP

plants = ["seattle", "san-diego"]
markets = ["new-york", "chicago", "topeka"]

capacity = [350, 600]
a = Dict(zip(plants, capacity))

demand = [325, 300, 275]
b = Dict(zip(markets, demand))

elasticity = [1.5, 1.2, 2.0]
esub = Dict(zip(markets, elasticity))

distance = [ 2.5 1.7 1.8 ;
             2.5 1.8 1.4  ]
d = Dict()
for i in 1:length(plants), j in 1:length(markets)
    d[plants[i], markets[j]] = distance[i,j]
end

f = 90

m = MCPModel()
@variable(m, w[i in plants] >= 0)
@variable(m, p[j in markets] >= 0)
@variable(m, x[i in plants, j in markets] >= 0)

@NLexpression(m, c[i in plants, j in markets], f * d[i,j] / 1000)

@operator(m, profit[i in plants, j in markets],    w[i] + c[i,j] - p[j])
@operator(m, supply[i in plants],                  a[i] - sum(x[i,j] for j in markets))
@operator(m, fxdemand[j in markets],               sum(x[i,j] for i in plants) - b[j])

@complementarity(m, profit, x)
@complementarity(m, supply, w)
@complementarity(m, fxdemand, p)

PATHSolver.options(convergence_tolerance=1e-8, output=:yes, time_limit=3600)


status = solveMCP(m)

@show getvalue(x)
@show getvalue(w)
@show getvalue(p)

@show status
@assert status == :Solved
@assert getvalue(x["seattle", "chicago"]) == 300.0
@assert getvalue(p["topeka"]) == 0.126
```

The result is
```julia
getvalue(x) = x: 2 dimensions:
[  seattle,:]
  [  seattle,new-york] = 49.99999533220467
  [  seattle, chicago] = 300.0
  [  seattle,  topeka] = 0.0
[san-diego,:]
  [san-diego,new-york] = 275.00000466779534
  [san-diego, chicago] = 0.0
  [san-diego,  topeka] = 275.0

getvalue(w) = w: 1 dimensions:
[  seattle] = 0.0
[san-diego] = 0.0

getvalue(p) = p: 1 dimensions:
[new-york] = 0.22499999999999992
[ chicago] = 0.15299999999999955
[  topeka] = 0.126

status = :Solved
```

## Status Symbols
```julia
status =
 [  :Solved,                          # 1 - solved
    :StationaryPointFound,            # 2 - stationary point found
    :MajorIterationLimit,             # 3 - major iteration limit
    :CumulativeMinorIterationLimit,   # 4 - cumulative minor iteration limit
    :TimeLimit,                       # 5 - time limit
    :UserInterrupt,                   # 6 - user interrupt
    :BoundError,                      # 7 - bound error (lb is not less than ub)
    :DomainError,                     # 8 - domain error (could not find a starting point)
    :InternalError                    # 9 - internal error
  ]
 ```


# Solution via NLsolve.jl

## Example 3

In the above Example 1, the only different part is to use `solver=:NLsolve`.
```julia
using Complementarity, JuMP
m = MCPModel()

lb = zeros(4)
ub = Inf*ones(4)
items = 1:4
@variable(m, lb[i] <= x[i in items] <= ub[i])

@operator(m, F1, 3*x[1]^2+2*x[1]*x[2]+2*x[2]^2+x[3]+3*x[4]-6)
@operator(m, F2, 2*x[1]^2+x[1]+x[2]^2+3*x[3]+2*x[4]-2)
@operator(m, F3, 3*x[1]^2+x[1]*x[2]+2*x[2]^2+2*x[3]+3*x[4]-1)
@operator(m, F4, x[1]^2+3*x[2]^2+2*x[3]+3*x[4]-3)

@complementarity(m, F1, x[1])
@complementarity(m, F2, x[2])
@complementarity(m, F3, x[3])
@complementarity(m, F4, x[4])

status = solveMCP(m, solver=:NLsolve, method=:trust_region)
@show status

z = getvalue(x)
@show z
```

The result should look like
```julia
status = Results of Nonlinear Solver Algorithm
 * Algorithm: Trust-region with dogleg and autoscaling
 * Starting Point: [0.0,0.0,0.0,0.0]
 * Zero: [1.2247448722522958,-7.67758278782679e-13,-9.147691538349567e-15,0.5000000002047756]
 * Inf-norm of residuals: 0.000000
 * Iterations: 132
 * Convergence: true
   * |x - x'| < 0.0e+00: false
   * |f(x)| < 1.0e-08: true
 * Function Calls (f): 133
 * Jacobian Calls (df/dx): 60

 z = x: 1 dimensions:
 [1] = 1.2247448722522958
 [2] = -7.67758278782679e-13
 [3] = -9.147691538349567e-15
 [4] = 0.5000000002047756
```

You can access the output of NLsolve by the following fieldnames
```julia
julia> fieldnames(status)
12-element Array{Symbol,1}:
 :method       
 :initial_x    
 :zero         
 :residual_norm
 :iterations   
 :x_converged  
 :xtol         
 :f_converged  
 :ftol         
 :trace        
 :f_calls      
 :g_calls      
```
For example:
```julia
julia> status.residual_norm
6.9373147226770016e-9

julia> status.x_converged
false

julia> status.f_converged
true
```
