using Complementarity, JuMP
using Base.Test

#########################################################################
m = MCPModel()

@defVar(m, x3 >= 0)
@defVar(m, x4 >= 0)
@defVar(m, x1 >= 0)
@defVar(m, x2 >= 0)

@defNLExpr(m, F2, x3-2x4 +2)
@defNLExpr(m, F3, x1-x2+2x3-2x4 -2)
@defNLExpr(m, F4, x1+2x2-2x3+4x4 -6)
@defNLExpr(m, F1, -x3-x4 +2)

complements(m, F4, x4)
complements(m, F1, x1)
complements(m, F3, x3)
complements(m, F2, x2)


PATHSolver.path_options(
                "convergence_tolerance 1e-8",
                "output no",
                "time_limit 3600"
                )
status = solveMCP(m)

z = [getValue(x1), getValue(x2), getValue(x3), getValue(x4)]
@show z


@test z == [2.8, 0.0, 0.8, 1.2]
#########################################################################


println("------------------------------------------------------------------")


#########################################################################
m = MCPModel()

@defVar(m, x3 >= 0)
@defVar(m, x4 >= 0)
@defVar(m, x1 >= 0)
@defVar(m, x2 >= 0)

@defNLExpr(m, F2, x3-2x4 +2)
@defNLExpr(m, F3, x1-x2+2x3-2x4 -2)
@defNLExpr(m, F4, x1+2x2-2x3+4x4 -6)
@defNLExpr(m, F1, -x3-x4 +2)

complements(m, F2, x2)
complements(m, F3, x3)
complements(m, F1, x1)
complements(m, F4, x4)

PATHSolver.path_options(
                "convergence_tolerance 1e-8",
                "output no",
                "time_limit 3600"
                )
status = solveMCP(m)

z = [getValue(x1), getValue(x2), getValue(x3), getValue(x4)]
@show z
@test z == [2.8, 0.0, 0.8, 1.2]
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

@defVar(m, lb[i] <= x[i in 1:4] <= ub[i])
@defNLExpr(m, F[i=1:4], sum{M[i,j]*x[j], j=1:4} + q[i])
complements(m, F, x)

PATHSolver.path_options(
                "convergence_tolerance 1e-8",
                "output no",
                "time_limit 3600"
                )

status = solveMCP(m)

z = getValue(x)
@test z[1] == 2.8
@test z[2] == 0.0
@test z[3] == 0.8
@test z[4] == 1.2
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

# @defVar(m, lb[i] <= x[i in items] <= ub[i])
@defVar(m, x[i in items] >= 0)
@defNLExpr(m, F[i in items], sum{M[i,j]*x[j], j in items} + q[i])
complements(m, F, x)

PATHSolver.path_options(
                "convergence_tolerance 1e-8",
                "output no",
                "time_limit 3600"
                )

status = solveMCP(m)

z = getValue(x)
@test z[1] == 2.8
@test z[2] == 0.0
@test z[3] == 0.8
@test z[4] == 1.2
#########################################################################
