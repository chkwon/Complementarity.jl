@testset "Updating bounds" begin
    m = MCPModel()
    @variable(m, x >= 0)
    @mapping(m, F, x+2)
    @complementarity(m, F, x)

    status = solveMCP(m)
    @test result_value(x) ≈ 0.

    JuMP.fix(x, 1., force=true)   
    status = solveMCP(m)
    @test result_value(x) ≈ 1.

    JuMP.unfix(x)
    JuMP.set_lower_bound(x, 0)
    status = solveMCP(m)
    @test result_value(x) ≈ 0.
end
