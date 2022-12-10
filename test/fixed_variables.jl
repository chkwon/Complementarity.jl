using Complementarity
using Test


@testset "Fixing Variables" begin



    A = (:a,:b,:c,:d)
    V = (a = 4, b = 3, c = 10, d = 2)

    #Step 1: Construct each mapping on the entire set A, and fix a variable. 
    m = MCPModel()
    @variable(m, x[A])
    @variable(m,y[A])
    @mapping(m, F[a = A], x[a]+V[a])
    @mapping(m, G[a = A], x[a]*y[a]+V[a])


    JuMP.fix(x[:b],2)

    @complementarity(m, F, x)
    @complementarity(m, G, y)

    status = solveMCP(m)


    @test all(result_value.(x) .≈ [-4,2,-10,-2])
    @test all(result_value.(y) .≈ [1,-1.5,1,1])

    @show result_value.(x)


    JuMP.fix(x[:b], 3)
    status = solveMCP(m)
    @test result_value.(x[:b]) ≈ 3


end

    #=
    
    #Step 2: Reconstruct the model, but restrict the domain of the equation with a fixed variable
    #A = (:a,:b,:c,:d)
    B = (:b)
    C = (:a,:c,:d)



    m = MCPModel()
    @variable(m, x[A])
    @variable(m,y[A])
    @mapping(m, F[a = C], x[a]+V[a])
    @mapping(m, G[a = A], x[a]*y[a]+V[a])    
    
    JuMP.fix(x[:b], 2., force=true)   
    
    @complementarity(m, F, x)
    @complementarity(m, G, y)
    
    status = solveMCP(m)

    @test all(result_value.(x) .≈ [-4,2,-10,-2])
    @test all(result_value.(y) .≈ [1,-1.5,1,1])


    #Step 3: Fix after the complementarity condition
    m = MCPModel()
    @variable(m, x[A])
    @variable(m,y[A])
    @mapping(m, F[a = C], x[a]+V[a])
    @mapping(m, G[a = A], x[a]*y[a]+V[a])    
    

    
    @complementarity(m, F, x)
    @complementarity(m, G, y)
    
    JuMP.fix(x[:b], 2., force=true)   

    status = solveMCP(m)

    @test all(result_value.(x) .≈ [-4,2,-10,-2])
    @test all(result_value.(y) .≈ [1,-1.5,1,1])



    #Step 4: Don't fix the variable, should raise error
    m = MCPModel()
    @variable(m, x[A])
    @variable(m,y[A])
    @mapping(m, F[a = C], x[a]+V[a])
    @mapping(m, G[a = A], x[a]*y[a]+V[a])    
    

    
    @complementarity(m, F, x)
    @complementarity(m, G, y)


    @test_throws AssertionError status = solveMCP(m)


    #Step 5: Test that the error is raised with NLsolve
    m = MCPModel()
    @variable(m, x[A])
    @variable(m,y[A])
    @mapping(m, F[a = C], x[a]+V[a])
    @mapping(m, G[a = A], x[a]*y[a]+V[a])    
    

    
    @complementarity(m, F, x)
    @complementarity(m, G, y)

    @test_throws AssertionError status = solveMCP(m)


end

=#