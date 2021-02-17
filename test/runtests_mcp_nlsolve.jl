using Complementarity
using Test

@info("-------[Testing Complementarity/NLsolve]------------------------------------------")

@testset "mcp test 0 with NLsolve example" begin

    lb = [0., 0., 0., 0.]
    ub = [Inf, Inf, Inf, Inf]
    x0 = [1.25, 0., 0., 0.5]

    m = nothing
    m = MCPModel()

    @variable(m, x[1:4] >= 0)
    @mapping(m, F1, 3*x[1]^2+2*x[1]*x[2]+2*x[2]^2+x[3]+3*x[4]-6)
    @mapping(m, F2, 2*x[1]^2+x[1]+x[2]^2+3*x[3]+2*x[4]-2)
    @mapping(m, F3, 3*x[1]^2+x[1]*x[2]+2*x[2]^2+2*x[3]+3*x[4]-1)
    @mapping(m, F4, x[1]^2+3*x[2]^2+2*x[3]+3*x[4]-3)
    @complementarity(m, [F1, F2, F3, F4], x)

    set_start_value.(x, x0)

    status = solveMCP(m, solver=:NLsolve)
    @show status

    z = result_value.(x)
    # Fz = [result_value(F1) result_value(F2) result_value(F3) result_value(F4)]

    @show z
    # @show Fz

    @test isapprox(z[1], 1.22474, atol=1e-4)
    @test isapprox(z[2], 0.0, atol=1e-4)
    @test isapprox(z[3], -1.378e-19, atol=1e-4)
    @test isapprox(z[4], 0.5, atol=1e-4)

    # @test isapprox(Fz[1], -1.26298e-9, atol=1e-4)
    # @test isapprox(Fz[2], 3.22474, atol=1e-4)
    # @test isapprox(Fz[3], 5.0, atol=1e-4)
    # @test isapprox(Fz[4], 3.62723e-11, atol=1e-4)



    set_start_value.(x, x0)
    status = solveMCP(m, solver=:PATH)
    @show status

    z = result_value.(x)
    # Fz = [result_value(F1) result_value(F2) result_value(F3) result_value(F4)]

    @show z
    # @show Fz

    @test isapprox(z[1], 1.22474, atol=1e-4)
    @test isapprox(z[2], 0.0, atol=1e-4)
    @test isapprox(z[3], -1.378e-19, atol=1e-4)
    @test isapprox(z[4], 0.5, atol=1e-4)

    # @test isapprox(Fz[1], -1.26298e-9, atol=1e-4)
    # @test isapprox(Fz[2], 3.22474, atol=1e-4)
    # @test isapprox(Fz[3], 5.0, atol=1e-4)
    # @test isapprox(Fz[4], 3.62723e-11, atol=1e-4)


    set_start_value.(x, x0)
    status = solveMCP(m, solver=:NLsolve)
    @show status

    z = result_value.(x)
    # Fz = [result_value(F1) result_value(F2) result_value(F3) result_value(F4)]

    @show z
    # @show Fz

    @test isapprox(z[1], 1.22474, atol=1e-4)
    @test isapprox(z[2], 0.0, atol=1e-4)
    @test isapprox(z[3], -1.378e-19, atol=1e-4)
    @test isapprox(z[4], 0.5, atol=1e-4)

    # @test isapprox(Fz[1], -1.26298e-9, atol=1e-4)
    # @test isapprox(Fz[2], 3.22474, atol=1e-4)
    # @test isapprox(Fz[3], 5.0, atol=1e-4)
    # @test isapprox(Fz[4], 3.62723e-11, atol=1e-4)


end

println("------------------------------------------------------------------")

@testset "mcp test 1 with NLsolve" begin
    m = nothing
    m = MCPModel()

    M = [0  0 -1 -1 ;
         0  0  1 -2 ;
         1 -1  2 -2 ;
         1  2 -2  4 ]

    q = [2; 2; -2; -6]

    lb = zeros(4)
    ub = Inf*ones(4)

    items = 1:4

    # @variable(m, lb[i] <= x[i in items] <= ub[i])
    @variable(m, x[i in items] >= 0)
    @mapping(m, F[i in items], sum(M[i,j]*x[j] for j in items) + q[i])
    @complementarity(m, F, x)

    status = solveMCP(m, solver=:NLsolve)
    @show status

    z = result_value.(x)
    # Fz = result_value(F)

    @show z
    # @show Fz

    @test isapprox(z[1], 2.8, atol=1e-4)
    @test isapprox(z[2], 0.0, atol=1e-4)
    @test isapprox(z[3], 0.8, atol=1e-4)
    @test isapprox(z[4], 1.2, atol=1e-4)
end

println("------------------------------------------------------------------")

@testset "mcp test 2 with NLsolve" begin

    m = nothing
    m = MCPModel()

    M = [0  0 -1 -1 ;
         0  0  1 -2 ;
         1 -1  2 -2 ;
         1  2 -2  4 ]

    q = [2; 2; -2; -6]

    lb = zeros(4)
    ub = Inf*ones(4)

    items = 1:4

    @variable(m, lb[i] <= x[i in items] <= ub[i])
    # @variable(m, x[i in items] >= 0)
    @mapping(m, F[i in items], sum(M[i,j]*x[j] for j in items) + q[i])
    @complementarity(m, F, x)

    status = solveMCP(m, solver=:NLsolve)
    @show status

    z = result_value.(x)
    # Fz = result_value(F) # currently produces an error

    @show z
    # @show Fz

    @test isapprox(z[1], 2.8, atol=1e-4)
    @test isapprox(z[2], 0.0, atol=1e-4)
    @test isapprox(z[3], 0.8, atol=1e-4)
    @test isapprox(z[4], 1.2, atol=1e-4)
end

println("------------------------------------------------------------------")

@testset "mcp test 3 with NLsolve" begin

    m = nothing
    m = MCPModel()

    M = [0  0 -1 -1 ;
         0  0  1 -2 ;
         1 -1  2 -2 ;
         1  2 -2  4 ]

    q = [2; 2; -2; -6]

    lb = zeros(4)
    ub = Inf*ones(4)

    @variable(m, lb[i] <= myvariablename[i in 1:4] <= ub[i])
    @mapping(m, myconst[i=1:4], sum(M[i,j]*myvariablename[j] for j in 1:4) + q[i])
    @complementarity(m, myconst, myvariablename)

    status = solveMCP(m, solver=:NLsolve)
    @show status

    z = result_value.(myvariablename)
    # Fz = result_value(myconst)

    @show z
    # @show Fz

    @test isapprox(z[1], 2.8, atol=1e-4)
    @test isapprox(z[2], 0.0, atol=1e-4)
    @test isapprox(z[3], 0.8, atol=1e-4)
    @test isapprox(z[4], 1.2, atol=1e-4)

end

println("------------------------------------------------------------------")

@testset "mcp test 4 with NLsolve" begin

    m = nothing
    m = MCPModel()

    lb = zeros(4)
    ub = Inf*ones(4)
    items = 1:4
    @variable(m, lb[i] <= x[i in items] <= ub[i])

    @mapping(m, F1, 3*x[1]^2 + 2*x[1]*x[2] + 2*x[2]^2 + x[3] + 3*x[4] -6)
    @mapping(m, F2, 2*x[1]^2 + x[1] + x[2]^2 + 3*x[3] + 2*x[4] -2)
    @mapping(m, F3, 3*x[1]^2 + x[1]*x[2] + 2*x[2]^2 + 2*x[3] + 3*x[4] -1)
    @mapping(m, F4, x[1]^2 + 3*x[2]^2 + 2*x[3] + 3*x[4] - 3)

    @complementarity(m, F1, x[1])
    @complementarity(m, F2, x[2])
    @complementarity(m, F3, x[3])
    @complementarity(m, F4, x[4])

    set_start_value(x[1], 1.25)
    set_start_value(x[2], 0.)
    set_start_value(x[3], 0.)
    set_start_value(x[4], 0.5)

    status = solveMCP(m, solver=:NLsolve)
    @show status

    z = result_value.(x)
    # Fz = [result_value(F1), result_value(F2), result_value(F3), result_value(F4)]

    @show z
    # @show Fz

    @test isapprox(z[1], 1.22474, atol=1e-4)
    @test isapprox(z[2], 0.0, atol=1e-4)
    @test isapprox(z[3], 0.0, atol=1e-4)
    @test isapprox(z[4], 0.5, atol=1e-4)




    set_start_value(x[1], 1.19)
    set_start_value(x[2], 0.2)
    set_start_value(x[3], 0.1)
    set_start_value(x[4], 0.49)

    status = solveMCP(m, solver=:PATH)
    @show status

    z = result_value.(x)
    # Fz = [result_value(F1), result_value(F2), result_value(F3), result_value(F4)]

    @show z
    # @show Fz

    @test isapprox(z[1], 1.22474, atol=1e-4)
    @test isapprox(z[2], 0.0, atol=1e-4)
    @test isapprox(z[3], 0.0, atol=1e-4)
    @test isapprox(z[4], 0.5, atol=1e-4)





    set_start_value(x[1], 1.40)
    set_start_value(x[2], 0.1)
    set_start_value(x[3], 0.2)
    set_start_value(x[4], 0.43)

    status = solveMCP(m, solver=:NLsolve)
    @show status

    z = result_value.(x)
    # Fz = [result_value(F1), result_value(F2), result_value(F3), result_value(F4)]

    @show z
    # @show Fz

    @test isapprox(z[1], 1.22474, atol=1e-4)
    @test isapprox(z[2], 0.0, atol=1e-4)
    @test isapprox(z[3], 0.0, atol=1e-4)
    @test isapprox(z[4], 0.5, atol=1e-4)
end
