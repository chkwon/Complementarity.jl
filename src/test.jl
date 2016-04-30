# using Complementarity, JuMP
using JuMP, Base.Test
import PATHSolver

include("mcp.jl")

m = MCPModel()

# @variable(m, 3 <= x3 <= 111)
# @variable(m, 4 <= x4 <= 222)
# @variable(m, 1 <= x1 <= 333)
# @variable(m, 2 <= x2 <= 444)

@variable(m, x3 >= 0)
@variable(m, x4 >= 0)
@variable(m, x1 >= 0)
@variable(m, x2 >= 0)

@NLexpression(m, F2, x3-2x4 +2)
@NLexpression(m, F3, x1-x2+2x3-2x4 -2)
@NLexpression(m, F4, x1+2x2-2x3+4x4 -6)
@NLexpression(m, F1, -x3-x4 +2)

# @variable(m, x3 >= 0)
# @variable(m, x4 >= 0)
# @variable(m, x1 >= 0)
# @variable(m, x2 >= 0)
#
# @function(m, F2, x3-2x4 +2)
# @function(m, F3, x1-x2+2x3-2x4 -2)
# @function(m, F4, x1+2x2-2x3+4x4 -6)
# @function(m, F1, -x3-x4 +2)

#
# @variable(m, x3 >= 30)
# @variable(m, x4 >= 40)
# @variable(m, x1 >= 10)
# @variable(m, x2 >= 20)
#
# @NLexpression(m, F2, x2*100)
# @NLexpression(m, F3, x3*100)
# @NLexpression(m, F4, x4*100)
# @NLexpression(m, F1, x1*100)


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

z = [getvalue(x1), getvalue(x2), getvalue(x3), getvalue(x4)]
@show z


@test z == [2.8, 0.0, 0.8, 1.2]
