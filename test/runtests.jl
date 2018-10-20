using Test

@testset "Complementarity Tests" begin

include("../examples/mcp/transmcp.jl")

include("hhoeschle_mcp1.jl")

include("runtests_mcp.jl")

include("runtests_lcp.jl")

include("runtests_mcp_nlsolve.jl")

include("runtests_mpec.jl")

end
