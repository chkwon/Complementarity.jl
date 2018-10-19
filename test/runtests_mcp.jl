using Complementarity, Test

@testset "mcp test 1 with PATHSolver" begin

    m = MCPModel()

    @variable(m, x3 >= 0)
    @variable(m, x4 >= 0)
    @variable(m, x1 >= 0)
    @variable(m, x2 >= 0)

    @mapping(m, F2, x3-2x4 +2)
    @mapping(m, F3, x1-x2+2x3-2x4 -2)
    @mapping(m, F4, x1+2x2-2x3+4x4 -6)
    @mapping(m, F1, -x3-x4 +2)

    @complementarity(m, F4, x4)
    @complementarity(m, F1, x1)
    @complementarity(m, F3, x3)
    @complementarity(m, F2, x2)


    PATHSolver.options(convergence_tolerance=1e-8, output=:yes, time_limit=3600)
    status = solveMCP(m, linear=true)

    z = [getvalue(x1), getvalue(x2), getvalue(x3), getvalue(x4)]
    Fz = [getvalue(F1), getvalue(F2), getvalue(F3), getvalue(F4)]
    @show z
    @show Fz

    @test isapprox(z, [2.8, 0.0, 0.8, 1.2])

    status = solveMCP(m)

    z = [getvalue(x1), getvalue(x2), getvalue(x3), getvalue(x4)]
    Fz = [getvalue(F1), getvalue(F2), getvalue(F3), getvalue(F4)]
    @show z
    @show Fz

    @test isapprox(z, [2.8, 0.0, 0.8, 1.2])

end


println("------------------------------------------------------------------")


@testset "mcp test 2 with PATHSolver" begin

    m = MCPModel()

    @variable(m, x3 >= 0)
    @variable(m, x4 >= 0)
    @variable(m, x1 >= 0)
    @variable(m, x2 >= 0)

    @mapping(m, F2, x3-2x4 +2)
    @mapping(m, F3, x1-x2+2x3-2x4 -2)
    @mapping(m, F4, x1+2x2-2x3+4x4 -6)
    @mapping(m, F1, -x3-x4 +2)

    @complementarity(m, F2, x2)
    @complementarity(m, F3, x3)
    @complementarity(m, F1, x1)
    @complementarity(m, F4, x4)

    PATHSolver.options(convergence_tolerance=1e-8, output=:yes, time_limit=3600)

    status = solveMCP(m)

    z = [getvalue(x1), getvalue(x2), getvalue(x3), getvalue(x4)]
    Fz = [getvalue(F1), getvalue(F2), getvalue(F3), getvalue(F4)]
    @show z
    @show Fz

    @test isapprox(z, [2.8, 0.0, 0.8, 1.2])
end


println("------------------------------------------------------------------")


@testset "mcp test 3 with PATHSolver" begin

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

    PATHSolver.options(convergence_tolerance=1e-8, output=:yes, time_limit=3600)


    status = solveMCP(m)

    z = getvalue(myvariablename)
    Fz = getvalue(myconst)

    @show z
    @show Fz

    @test isapprox(z[1], 2.8)
    @test isapprox(z[2], 0.0)
    @test isapprox(z[3], 0.8)
    @test isapprox(z[4], 1.2)
end

println("------------------------------------------------------------------")


@testset "mcp test 4 with PATHSolver" begin

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

    PATHSolver.options(convergence_tolerance=1e-8, output=:no, time_limit=3600)


    status = solveMCP(m)

    z = getvalue(x)
    # Fz = getvalue(F) # currently produces an error

    @show z
    # @show Fz

    @test isapprox(z[1], 2.8)
    @test isapprox(z[2], 0.0)
    @test isapprox(z[3], 0.8)
    @test isapprox(z[4], 1.2)
end


println("------------------------------------------------------------------")

@testset "mcp test 5 with PATHSolver" begin

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

    setvalue(x[1], 1.25)
    setvalue(x[2], 0.)
    setvalue(x[3], 0.)
    setvalue(x[4], 0.5)

    status = solveMCP(m)
    @show status

    z = getvalue(x)
    Fz = [getvalue(F1), getvalue(F2), getvalue(F3), getvalue(F4)]

    @show z
    @show Fz

    @test isapprox(z[1], 1.22474, atol=1e-4)
    @test isapprox(z[2], 0.0, atol=1e-4)
    @test isapprox(z[3], 0.0, atol=1e-4)
    @test isapprox(z[4], 0.5, atol=1e-4)


    setvalue(x[1], 1.19)
    setvalue(x[2], 0.2)
    setvalue(x[3], 0.1)
    setvalue(x[4], 0.49)

    status = solveMCP(m)
    @show status

    z = getvalue(x)
    Fz = [getvalue(F1), getvalue(F2), getvalue(F3), getvalue(F4)]

    @show z
    @show Fz

    @test isapprox(z[1], 1.22474, atol=1e-4)
    @test isapprox(z[2], 0.0, atol=1e-4)
    @test isapprox(z[3], 0.0, atol=1e-4)
    @test isapprox(z[4], 0.5, atol=1e-4)
end


println("------------------------------------------------------------------")
