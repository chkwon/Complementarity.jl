# using Complementarity, JuMP
using JuMP, Base.Test
import PATHSolver

include("mcp.jl")

m = MCPModel()

# @defVar(m, 3 <= x3 <= 111)
# @defVar(m, 4 <= x4 <= 222)
# @defVar(m, 1 <= x1 <= 333)
# @defVar(m, 2 <= x2 <= 444)

@defVar(m, x3 >= 0)
@defVar(m, x4 >= 0)
@defVar(m, x1 >= 0)
@defVar(m, x2 >= 0)

@defNLExpr(m, F2, x3-2x4 +2)
@defNLExpr(m, F3, x1-x2+2x3-2x4 -2)
@defNLExpr(m, F4, x1+2x2-2x3+4x4 -6)
@defNLExpr(m, F1, -x3-x4 +2)


#
# @defVar(m, x3 >= 30)
# @defVar(m, x4 >= 40)
# @defVar(m, x1 >= 10)
# @defVar(m, x2 >= 20)
#
# @defNLExpr(m, F2, x2*100)
# @defNLExpr(m, F3, x3*100)
# @defNLExpr(m, F4, x4*100)
# @defNLExpr(m, F1, x1*100)


complements(m, F2, x2)
complements(m, F3, x3)
complements(m, F1, x1)
complements(m, F4, x4)





PATHSolver.path_options(
                "convergence_tolerance 1e-8",
                "output no",
                "time_limit 3600"
                )
solveMCP(m)

z = [getValue(x1), getValue(x2), getValue(x3), getValue(x4)]
@show z


@test z == [2.8, 0.0, 0.8, 1.2]
