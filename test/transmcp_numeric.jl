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

using Complementarity
using Test

@testset "transmcp_numeric.jl 1" begin

    # plants = ["seattle", "san-diego"]
    plants = 1:2
    # markets = ["new-york", "chicago", "topeka"]
    markets = 1:3

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
    @test status == :Solved
    @test isapprox(result_value(x[1,2]), 300.0)
    @test isapprox(result_value(p[3]),  0.126)

end



@testset "transmcp_numeric.jl 2" begin

    # plants = ["seattle", "san-diego"]
    plants = 1:2
    # markets = ["new-york", "chicago", "topeka"]
    markets = 1:3

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
    @variable(m, w[1:2] >= 0)
    @variable(m, p[1:3] >= 0)
    @variable(m, x[1:2, 1:3] >= 0)

    @NLexpression(m, c[i in 1:2, j in 1:3], f * d[i,j] / 1000)
    @mapping(m, profit[i in 1:2, j in 1:3], w[i] + c[i,j] - p[j])
    @mapping(m, supply[i in 1:2], a[i] - sum(x[i,j] for j in markets))
    @mapping(m, fxdemand[j in 1:3], sum(x[i,j] for i in plants) - b[j])

    @complementarity(m, profit, x)
    @complementarity(m, supply, w)
    @complementarity(m, fxdemand, p)

    status = solveMCP(m; convergence_tolerance=1e-8, output="yes", time_limit=3600)


    @show result_value.(x)
    @show result_value.(w)
    @show result_value.(p)

    @show status
    @test status == :Solved
    @test isapprox(result_value(x[1,2]), 300.0)
    @test isapprox(result_value(p[3]),  0.126)

end
