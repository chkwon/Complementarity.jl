# # dempe.mod	QQR2-MN-4-3
# # Original AMPL coding by Sven Leyffer
#
# # An MPEC from S. Dempe, "A necessary and sufficient optimality
# # condition for bilevel programming problems", Optimization 25,
# # pp. 341-354, 1992. From book Math. Progr. with Equil. Constr,
# # by Luo, Pang & Ralph, CUP, 1997, p. 354.
#
# # Number of variables:   2 + 1 multipliers
# # Number of constraints: 2
# # Nonlinear complementarity constraints
#
# var x;
# var z;
# var w >= 0;
#
# minimize f: (x - 3.5)^2 + (z + 4)^2;
#
# subject to
#     con1:  z - 3 + 2*z*w = 0;
#     con2:  0 >= z^2 - x     complements  w >= 0;
#
# data;
#
# # starting point
# let x := 1;
# let z := 1;
# let w := 1;
#
#
# let x :=   0.183193 ;
# let z := 0.428106;
# let w := 3.00379;


using JuMP, Complementarity
using Ipopt
using Test

@testset "dempe.jl" begin

    dempe = Model(Ipopt.Optimizer)


    @variable(dempe, x)
    @variable(dempe, z)
    @variable(dempe, w>=0)

    @NLobjective(dempe, Min, (x - 3.5)^2 + (z + 4)^2)

    @NLconstraint(dempe, z - 3 + 2*z*w == 0)

    @complements(dempe, 0 >= z^2 - x,  w >= 0, :simple)

    # Initial solutions to help reaching the optimality
    set_start_value(x, 0)
    set_start_value(z, 0)
    set_start_value(w, 1e7)

    # (xx, zz, ww) = (-3.999779029254958e-10, 9.852984581857451e-8, 1.5223813031201378e7)


    JuMP.optimize!(dempe)

    @show JuMP.objective_value(dempe)
    @test isapprox(JuMP.objective_value(dempe), 28.25, atol=1e-4)

    xx = JuMP.value(x)
    zz = JuMP.value(z)
    ww = JuMP.value(w)
    @test isapprox(zz - 3 + 2*zz*ww, 0, atol=1e-4)
end
