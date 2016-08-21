# using Complementarity, JuMP
using JuMP, Base.Test
import PATHSolver, NLsolve

include("mcp.jl")


#########################################################################
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
@NLexpression(m, F[i in items], sum{M[i,j]*x[j], j in items} + q[i])
complements(m, F, x)

status = solveMCP(m, solver=:NLsolve)
@show status

z = getvalue(x)
@show z
# @test isapprox(z[1], 2.8)
# @test isapprox(z[2], 0.0)
# @test isapprox(z[3], 0.8)
# @test isapprox(z[4], 1.2)
#########################################################################


#########################################################################
m = MCPModel()

items = 1:4
@variable(m, x[i in items] >= 0)

@NLexpression(m, F1, 3*x[1]^2+2*x[1]*x[2]+2*x[2]^2+x[3]+3*x[4]-6)
@NLexpression(m, F2, 2*x[1]^2+x[1]+x[2]^2+3*x[3]+2*x[4]-2)
@NLexpression(m, F3, 3*x[1]^2+x[1]*x[2]+2*x[2]^2+2*x[3]+3*x[4]-1)
@NLexpression(m, F4, x[1]^2+3*x[2]^2+2*x[3]+3*x[4]-3)

complements(m, F1, x[1])
complements(m, F2, x[2])
complements(m, F3, x[3])
complements(m, F4, x[4])

status = solveMCP(m, solver=:NLsolve, method=:trust_region)
@show status

z = getvalue(x)
@show z

@test isapprox(z[1], 1.2247, atol=1e-4)
@test isapprox(z[2], 0.0, atol=1e-4)
@test isapprox(z[3], 0.0, atol=1e-4)
@test isapprox(z[4], 0.5, atol=1e-4)

#########################################################################
