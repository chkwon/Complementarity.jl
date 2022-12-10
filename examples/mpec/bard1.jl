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


# @testset "bard1.jl" begin
    bard1 = Model(Ipopt.Optimizer)

    @variable(bard1, x>=0)
    @variable(bard1, y>=0)
    @variable(bard1, l[1:3])

    @NLobjective(bard1, Min, (x - 5)^2 + (2*y + 1)^2)

    @NLconstraint(bard1, 2*(y-1) - 1.5*x + l[1] - l[2]*0.5 + l[3] == 0)

    @complements(bard1, 0 <= 3*x - y - 3,        l[1] >= 0, :smooth)
    @complements(bard1, 0 <= - x + 0.5*y + 4,    l[2] >= 0, :smooth)
    @complements(bard1, 0 <= - x - y + 7,        l[3] >= 0, :simple)

    JuMP.optimize!(bard1)

    @show JuMP.objective_value(bard1)
    @test isapprox(JuMP.objective_value(bard1), 17.0000, atol=1e-4)
# end