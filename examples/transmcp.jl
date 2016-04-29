# transmcp.jl
# translated by Changhyun Kwon from transmcp.gms
# http://www.gams.com/modlib/libhtml/transmcp.htm
#
# Transportation model as equilibrium problem (TRANSMCP,SEQ=126)
#
#    Dantzig's original transportation model (TRNSPORT) is
#    reformulated as a linear complementarity problem.  We first
#    solve the model with fixed demand and supply quantities, and
#    then we incorporate price-responsiveness on both sides of the
#    market.
#
#
# Dantzig, G B, Chapter 3.3. In Linear Programming and Extensions.
# Princeton University Press, Princeton, New Jersey, 1963.

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

using Complementarity, JuMP

m = MCPModel()
@defVar(m, w[i in plants] >= 0)
@defVar(m, p[j in markets] >= 0)
@defVar(m, x[i in plants, j in markets] >= 0)

@defNLExpr(m, c[i in plants, j in markets], f * d[i,j] / 1000)

@defNLExpr(m, profit[i in plants, j in markets], w[i] + c[i,j] - p[j])
@defNLExpr(m, supply[i in plants], a[i] - sum{x[i,j], j in markets})
@defNLExpr(m, fxdemand[j in markets], sum{x[i,j], i in plants} - b[j])

correspond(m, profit, x)
correspond(m, supply, w)
correspond(m, fxdemand, p)

PATHSolver.path_options(
                "convergence_tolerance 1e-2",
                "output no",
                "time_limit 3600"
                )

solveMCP(m)

@show getValue(x)
@show getValue(w)
@show getValue(p)
