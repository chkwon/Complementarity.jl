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
When `F(x)` is a linear mapping such as `F(x) = M x + q` with matrix `M` and vector `q`, then it is a Linear Complementarity Problem (LCP). All these problems are solved by the [PATH Solver](http://pages.cs.wisc.edu/%7Eferris/path.html) which is wrapped by the [PATHSolver.jl](https://github.com/chkwon/PATHSolver.jl) package.

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
@mapping(m, F[i in items], sum(M[i,j]*x[j] for j in items) + q[i])
@complementarity(m, F, x)


status = solveMCP(m, solver=:PATH, convergence_tolerance=1e-8, output="yes", time_limit=3600)

z = result_value.(x)
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
@mapping(m, F[i in items], sum(M[i,j]*x[j] for j in items) + q[i])
```
This is to define expressions for `F` in MCP. This is merely an alias of `JuMP.@NLexpression`.

```julia
@complementarity(m, F, x)
```
This macro matches each element of `F` and the complementing element of `x`.


```julia
status = solveMCP(m, solver=:PATH, convergence_tolerance=1e-8, output="yes", time_limit=3600)
```
This solves the MCP and stores the solution inside `m`, which can be accessed by `result_value(x)`.
Keyword arguments are options of the PATH Solver. See the [list of options](http://www.cs.wisc.edu/~ferris/path/options.pdf).




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

@mapping(m, profit[i in plants, j in markets],    w[i] + c[i,j] - p[j])
@mapping(m, supply[i in plants],                  a[i] - sum(x[i,j] for j in markets))
@mapping(m, fxdemand[j in markets],               sum(x[i,j] for i in plants) - b[j])

@complementarity(m, profit, x)
@complementarity(m, supply, w)
@complementarity(m, fxdemand, p)

status = solveMCP(m; convergence_tolerance=1e-8, output="yes", time_limit=3600)

@show result_value.(x)
@show result_value.(w)
@show result_value.(p)

@show status
@assert status == :Solved
@assert result_value(x["seattle", "chicago"]) == 300.0
@assert result_value(p["topeka"]) == 0.126
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

## Warmstart

When you need warmstart, you can do either:

```julia
x_start = Dict(
    ("seattle","new-york")=>50,
    ("seattle","chicago")=>200,
    ("seattle","topeka")=>0,
    ("san-diego","new-york")=>275,
    ("san-diego","chicago")=>0,
    ("san-diego","topeka")=>275,
    )
@variable(m, x[i in plants, j in markets] >= 0, start=x_start[(i,j)])
```

or

```julia
@variable(m, x[i in plants, j in markets] >= 0)
setvalue(x["seattle", "chicago"], 200)
```



## Termination Status
```julia
    :Solved,                          # 1 - The problem was solved
    :StationaryPointFound,            # 2 - A stationary point was found
    :MajorIterationLimit,             # 3 - Major iteration limit met
    :CumulativeMinorIterationLimit,   # 4 - Cumulative minor iterlim met
    :TimeLimit,                       # 5 - Ran out of time
    :UserInterrupt,                   # 6 - Control-C, typically
    :BoundError,                      # 7 - Problem has a bound error (lb is not less than ub)
    :DomainError,                     # 8 - Could not find starting point
    :Infeasible,                      # 9 - Problem has no solution  
    :Error,                           #10 - An error occurred within the code
    :LicenseError,                    #11 - License could not be found
    :OK                               #12 - OK
```

# Solution via NLsolve.jl

## Example 3

We can specify `NLsolve` as the solver by providing `solver=:NLsolve` to `solveMCP()`.
```julia
using Complementarity

m = MCPModel()

lb = zeros(4)
ub = Inf*ones(4)
items = 1:4
@variable(m, lb[i] <= x[i in items] <= ub[i])

@mapping(m, F1, 3*x[1]^2 + 2*x[1]*x[2] + 2*x[2]^2 + x[3] + 3*x[4] -6)
@mapping(m, F2, 2*x[1]^2 + x[1] + x[2]^2 + 3*x[3] + 2*x[4] -2)
@mapping(m, F3, 3*x[1]^2 + x[1]*x[2] + 2*x[2]^2 + 2*x[3] + 3*x[4] -1)
@mapping(m, F4, x[1]^2 + 3*x[2]^2 + 2*x[3] + 3*x[4] - 3)

@complementarity(m, F1, x[1])
@complementarity(m, F2, x[2])
@complementarity(m, F3, x[3])
@complementarity(m, F4, x[4])

setvalue(x[1], 1.25)
setvalue(x[2], 0.)
setvalue(x[3], 0.)
setvalue(x[4], 0.5)

status = solveMCP(m, solver=:NLsolve)
@show status

z = result_value(x)

@show z
@show Fz
```

Note that initial values are provided using `setvalue()`.

The result should look like
```julia
status = Results of Nonlinear Solver Algorithm
 * Algorithm: Trust-region with dogleg and autoscaling
 * Starting Point: [1.25,0.0,0.0,0.5]
 * Zero: [1.22474,0.0,-2.02379e-19,0.5]
 * Inf-norm of residuals: 0.000000
 * Iterations: 3
 * Convergence: true
   * |x - x'| < 0.0e+00: false
   * |f(x)| < 1.0e-08: true
 * Function Calls (f): 4
 * Jacobian Calls (df/dx): 4

z = x: 1 dimensions:
[1] = 1.2247448711263813
[2] = 0.0
[3] = -2.0237901522246342e-19
[4] = 0.5000000002286319

Fz = [-1.26298e-9,3.22474,5.0,3.62723e-11]
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
1.2629755019588629e-9

julia> status.x_converged
false

julia> status.f_converged
true
```
