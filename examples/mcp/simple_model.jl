using Complementarity

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

#m = JuMP.Model()

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


mcp_data = Complementarity.get_MCP_data(m)
p = Complementarity.sortperm_MCP_data(mcp_data)
n = length(mcp_data)


#JuMP.@NLconstraint(m,[i=1:n],mcp_data[p[i]].F == 0)