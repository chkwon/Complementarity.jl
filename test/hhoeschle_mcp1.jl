#########################################################################
# https://github.com/chkwon/PATHSolver.jl/issues/16
#########################################################################

using Complementarity, Printf

@testset "hhoeschle mcp 1" begin

    model = MCPModel()
    scenarios = [@sprintf "s%03d" s for s=1:3]
    INDEX = Dict("s001" => 1., "s002" => 0., "s003" => 0.5)
    COST = Dict("s001" => 10., "s002" => 2., "s003" => 5)
    COST_INDEX = Dict("s001" => 2., "s002" => 1., "s003" => 3.)

    @variable(model, xname[s in scenarios] >= 0)
    @variable(model, ynoname[s in scenarios; INDEX[s] >0] >= 0)
    @variable(model, lambda1 >= 0)
    @variable(model, lambda2 >= 0)
    @NLexpression(model, lambda, lambda1 - lambda2)

    @mapping(model, market_lambda1,
        + sum(xname[s] for s in scenarios)
        + sum(ynoname[s] for s in scenarios if INDEX[s] > 0)
        - 10)
    @complementarity(model, market_lambda1, lambda1)

    @mapping(model, market_lambda2,
        - sum(xname[s] for s in scenarios)
        - sum(ynoname[s] for s in scenarios if INDEX[s] > 0)
        + 10)
    @complementarity(model, market_lambda2, lambda2)

    @mapping(model, KKT_xname[s in scenarios],
        - lambda
        + COST[s])
    @complementarity(model, KKT_xname, xname)


    @mapping(model, KKT_ynoname[s in scenarios; INDEX[s] >0],
        - lambda
        + COST_INDEX[s])
    # @complementarity(model, KKT_ynoname, ynoname)
    for s in scenarios
        if INDEX[s] > 0
            @complementarity(model, KKT_ynoname[s], ynoname[s])
        end
    end



    stat = solveMCP(
        model;         
        convergence_tolerance=1e-6,
        output="yes",
        time_limit=3600*12,
        # time_limit=0.1,
        minor_iteration_limit=1000,
        cumulative_iteration_limit=1000000000
    )

    stat = solveMCP(
        model;         
        convergence_tolerance=1e-6,
        output="yes",
        time_limit=3600*12,
        # time_limit=0.1,
        minor_iteration_limit=1000,
        cumulative_iteration_limit=1000000000
    ) # solve twice to see if model is reusable

    stat = solveMCP(
        model;         
        convergence_tolerance=1e-6,
        output="yes",
        time_limit=3600*12,
        # time_limit=0.1,
        minor_iteration_limit=1000,
        cumulative_iteration_limit=1000000000
    ) # solve three times to see if model is reusable
end
