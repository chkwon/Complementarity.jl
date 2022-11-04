using Test
using Complementarity

@testset "Complementarity Tests" begin

    include("../examples/mcp/transmcp.jl")
    include("transmcp_numeric.jl")

    @testset "Array{Variable,4} Test" begin
        m = MCPModel()
        @variable(m, x[1:2, 1:3, 1:4, 1:5] >=0 )
        @mapping(m, F[i=1:2, j=1:3, k=1:4, l=1:5], x[i,j,k,l])
        # for i=1:2, j=1:3, k=1:4, l=1:5
        #     @complementarity(m, F[i,j,k,l], x[i,j,k,l])
        # end
        @complementarity(m, F, x)        
        status = solveMCP(m)
        @test status == :Solved
    end

    include("hhoeschle_mcp1.jl")

    include("runtests_mcp.jl")

    include("runtests_mcp_nlsolve.jl")

    include("runtests_mpec.jl")

    include("test_bounds.jl")

    include("fixed_variables.jl")
end
