using Complementarity, JuMP
using Base.Test

#########################################################################
m = MCPModel()

@defVar(m, x1)
@defVar(m, x2)
@defVar(m, x3)
@defVar(m, x4)

@defNLExpr(m, F1, -x3-x4 +2)
@defNLExpr(m, F2, x3-2x4 +2)
@defNLExpr(m, F3, x1-x2+2x3-2x4 -2)
@defNLExpr(m, F4, x1+2x2-2x3+4x4 -6)

correspond(m, 0, x1, Inf, F1)
correspond(m, 0, x2, Inf, F2)
correspond(m, 0, x3, Inf, F3)
correspond(m, 0, x4, Inf, F4)

PATHSolver.path_options(
                "convergence_tolerance 1e-2",
                "output no",
                "time_limit 3600"
                )
solveMCP(m)

z = [getValue(x1), getValue(x2), getValue(x3), getValue(x4)]
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

@defVar(m, x[1:4])
@defNLExpr(m, F[i=1:4], sum{M[i,j]*x[j], j=1:4} + q[i])
correspond(m, lb, x, ub, F)

PATHSolver.path_options(
                "convergence_tolerance 1e-2",
                "output yes",
                "time_limit 3600"
                )

solveMCP(m)

z = [getValue(x1), getValue(x2), getValue(x3), getValue(x4)]
@test z == [2.8, 0.0, 0.8, 1.2]
#########################################################################
