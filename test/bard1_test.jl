# # bard1.mod	QQR2-MN-8-5
# # Original AMPL coding by Sven Leyffer
#
# # An MPEC from J.F. Bard, Convex two-level optimization,
# # Mathematical Programming 40(1), 15-27, 1988.
#
# # Number of variables:   2 + 3 slack + 3 multipliers
# # Number of constraints: 4
#
# var x >= 0;
# var y >= 0;
#
# # ... multipliers
# var l{1..3};
#
# minimize f:(x - 5)^2 + (2*y + 1)^2;
#
# subject to
#
#    KKT:    2*(y-1) - 1.5*x + l[1] - l[2]*0.5 + l[3] = 0;
#
#    lin_1:  0 <= 3*x - y - 3        complements l[1] >= 0;
#    lin_2:  0 <= - x + 0.5*y + 4    complements l[2] >= 0;
#    lin_3:  0 <= - x - y + 7        complements l[3] >= 0;


using JuMP, Complementarity
using Ipopt
using Test


@testset "bard1_test.jl" begin

    m = Model(Ipopt.Optimizer)

    @variable(m, x>=0)
    @variable(m, y>=0)
    @variable(m, l[1:3])

    @NLobjective(m, Min, (x - 5)^2 + (2*y + 1)^2)

    # @NLconstraint(m, 2*(y-1) - 1.5*x + l[1] - l[2]*0.5 + l[3] == 0)
    @NLconstraint(m, 2*(y-1) - 1.5*x + l[1] - l[2]*0.5 + l[3] == 0)

    # Original Form
    # @complements(m, 0 <= 3*x - y - 3,        l[1] >= 0)
    # @complements(m, 0 <= - x + 0.5*y + 4,    l[2] >= 0)
    # @complements(m, 0 <= - x - y + 7,        l[3] >= 0)
    @complements(m, 0 >= - 3*x + y + 3,        l[1] >= 0, :smooth)
    @complements(m, 0 <= - x + 0.5*y + 4,      l[2] <= 0, :simple)
    @complements(m, 0 >= + x + y - 7,          l[3] <= 0, :smooth)

    # just testing
    @macroexpand(@complements(m, 3*x - y - 3,    0 <=  l[1] <= 5))
    @macroexpand(@complements(m, 3*x - y - 3,    5 >=  l[1] >= 0))
    
    @test_throws LoadError @macroexpand(@complements(m,  0 <= 3*x - y - 3 <=10,  l[1]))
    @test_throws LoadError @macroexpand(@complements(m, 10 >= 3*x - y - 3 >= 5,  l[1]))

    JuMP.optimize!(m)
    @show JuMP.value.(x)
    @show JuMP.value.(y)
    @show JuMP.value.(l)
    @show JuMP.objective_value(m)
    @test isapprox(JuMP.objective_value(m), 17.0000, atol=1e-4)
end
