using Complementarity, JuMP
using Base.Test

include("../examples/mcp/transmcp.jl")


#########################################################################
m = MCPModel()

@variable(m, x3 >= 0)
@variable(m, x4 >= 0)
@variable(m, x1 >= 0)
@variable(m, x2 >= 0)

@operator(m, F2, x3-2x4 +2)
@operator(m, F3, x1-x2+2x3-2x4 -2)
@operator(m, F4, x1+2x2-2x3+4x4 -6)
@operator(m, F1, -x3-x4 +2)

@complementarity(m, F4, x4)
@complementarity(m, F1, x1)
@complementarity(m, F3, x3)
@complementarity(m, F2, x2)


PATHSolver.path_options(
                "convergence_tolerance 1e-8",
                "output yes",
                "time_limit 3600"
                )
status = solveMCP(m)

z = [getvalue(x1), getvalue(x2), getvalue(x3), getvalue(x4)]
@show z


@test isapprox(z, [2.8, 0.0, 0.8, 1.2])


#########################################################################


println("------------------------------------------------------------------")


#########################################################################
m = MCPModel()

@variable(m, x3 >= 0)
@variable(m, x4 >= 0)
@variable(m, x1 >= 0)
@variable(m, x2 >= 0)

@operator(m, F2, x3-2x4 +2)
@operator(m, F3, x1-x2+2x3-2x4 -2)
@operator(m, F4, x1+2x2-2x3+4x4 -6)
@operator(m, F1, -x3-x4 +2)

@complementarity(m, F2, x2)
@complementarity(m, F3, x3)
@complementarity(m, F1, x1)
@complementarity(m, F4, x4)

PATHSolver.path_options(
                "convergence_tolerance 1e-8",
                "output yes",
                "time_limit 3600"
                )
status = solveMCP(m)

z = [getvalue(x1), getvalue(x2), getvalue(x3), getvalue(x4)]
@show z
@test isapprox(z, [2.8, 0.0, 0.8, 1.2])
#########################################################################


println("------------------------------------------------------------------")


#########################################################################
m = MCPModel()

M = [0  0 -1 -1 ;
     0  0  1 -2 ;
     1 -1  2 -2 ;
     1  2 -2  4 ]

q = [2; 2; -2; -6]

lb = zeros(4)
ub = Inf*ones(4)

@variable(m, lb[i] <= myvariablename[i in 1:4] <= ub[i])
@operator(m, myconst[i=1:4], sum(M[i,j]*myvariablename[j] for j in 1:4) + q[i])
@complementarity(m, myconst, myvariablename)

PATHSolver.path_options(
                "convergence_tolerance 1e-8",
                "output yes",
                "time_limit 3600"
                )

status = solveMCP(m)

z = getvalue(myvariablename)
@test isapprox(z[1], 2.8)
@test isapprox(z[2], 0.0)
@test isapprox(z[3], 0.8)
@test isapprox(z[4], 1.2)
#########################################################################

println("------------------------------------------------------------------")


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

@variable(m, lb[i] <= x[i in items] <= ub[i])
# @variable(m, x[i in items] >= 0)
@operator(m, F[i in items], sum(M[i,j]*x[j] for j in items) + q[i])
@complementarity(m, F, x)

PATHSolver.path_options(
                "convergence_tolerance 1e-8",
                "output no",
                "time_limit 3600"
                )

status = solveMCP(m)

z = getvalue(x)
@test isapprox(z[1], 2.8)
@test isapprox(z[2], 0.0)
@test isapprox(z[3], 0.8)
@test isapprox(z[4], 1.2)
#########################################################################


println("------------------------------------------------------------------")

println("-------[Testing NLSolve]------------------------------------------")


#########################################################################
m = nothing
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
@operator(m, F[i in items], sum(M[i,j]*x[j] for j in items) + q[i])
@complementarity(m, F, x)

status = solveMCP(m, solver=:NLsolve)
@show status

z = getvalue(x)
@show z
# @test isapprox(z[1], 2.8)
# @test isapprox(z[2], 0.0)
# @test isapprox(z[3], 0.8)
# @test isapprox(z[4], 1.2)
#########################################################################

println("------------------------------------------------------------------")

#########################################################################
m = nothing
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

# @test isapprox(z[1], 1.2247, atol=1e-4)
# @test isapprox(z[2], 0.0, atol=1e-4)
# @test isapprox(z[3], 0.0, atol=1e-4)
# @test isapprox(z[4], 0.5, atol=1e-4)

#########################################################################
