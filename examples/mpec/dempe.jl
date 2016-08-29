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

m = Model(solver=IpoptSolver())


@variable(m, x)
@variable(m, z)
@variable(m, w>=0)

@NLobjective(m, Min, (x - 3.5)^2 + (z + 4)^2)

@NLconstraint(m, z - 3 + 2*z*w == 0)

@complements(m, 0 >= z^2 - x,  w >= 0, simple)

# setvalue(x, 50)
# setvalue(z, 50)
# setvalue(w, 1e5)

solve(m)

@show getobjectivevalue(m)
@assert isapprox(getobjectivevalue(m), 28.25, atol=1e-4)
