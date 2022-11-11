# # data file for gnash1.mod
# # Original AMPL coding by Sven Leyffer, University of Dundee
#
# # An MPEC from F. Facchinei, H. Jiang and L. Qi, A smoothing method for
# # mathematical programs with equilibrium constraints, Universita di Roma
# # Technical report, 03.96. Problem number 8
#
# param	:	c	K	b	:=
# 	1	10	5	1.2
# 	2	 8	5	1.1
# 	3	 6	5	1.0
# 	4	 4	5	0.9
# 	5	 2	5	0.8	;
#
# param L := 20;
# param g := 1.7;
#
# let x := 75;


c = [10, 8, 6, 4, 2]
K = [5, 5, 5, 5, 5]
b = [1.2, 1.1, 1.0, 0.9, 0.8]
L = 20
g = 1.7
# x = 75


# # gnash1m.mod	QQR2-MN-10-8
# # Original AMPL coding by Sven Leyffer, University of Dundee
#
# # Formulation with mixed complementarity of ...
# #
# # An MPEC from F. Facchinei, H. Jiang and L. Qi, A smoothing method for
# # mathematical programs with equilibrium constraints, Universita di Roma
# # Technical report, 03.96. Problem number 8
# #
# # Arises from Gournot Nash equilibrium , 10 instances are available
# # (see gnash1i.dat, i=0,1,...,9)
#
# # Number of variables:   10
# # Number of constraints:  8
#
# # ... parameters for each firm/company
# param c{1..5};			# c_i
# param K{1..5};			# K_i
# param b{1..5};			# \beta_i
#
# # ... parameters for each problem instance
# param L;			# L
# param g;			# \gamma
#
# # ... computed constants
# param gg := 5000^(1/g);
#
# var x >= 0, <= L;
# var y{1..4};
# var l{1..4};   		        # Multipliers
# var Q = x+y[1]+y[2]+y[3]+y[4];	# defined variable Q
#
# minimize f: c[1]*x + b[1]/(b[1]+1)*K[1]^(-1/b[1])*x^((1+b[1])/b[1])
# 		- x*( gg*Q^(-1/g) );
#
# subject to
#
#    F1: 0 = ( c[2] + K[2]^(-1/b[2])*y[1] ) - ( gg*Q^(-1/g) )
# 				- y[1]*( -1/g*gg*Q^(-1-1/g) ) - l[1];
#    F2: 0 = ( c[3] + K[3]^(-1/b[3])*y[2] ) - ( gg*Q^(-1/g) )
# 				- y[2]*( -1/g*gg*Q^(-1-1/g) ) - l[2];
#    F3: 0 = ( c[4] + K[4]^(-1/b[4])*y[3] ) - ( gg*Q^(-1/g) )
# 				- y[3]*( -1/g*gg*Q^(-1-1/g) ) - l[3];
#    F4: 0 = ( c[5] + K[5]^(-1/b[5])*y[4] ) - ( gg*Q^(-1/g) )
# 				- y[4]*( -1/g*gg*Q^(-1-1/g) ) - l[4];
#
#    g1: 0 <= y[1] <= L  complements   l[1];
#    g3: 0 <= y[2] <= L  complements   l[2];
#    g5: 0 <= y[3] <= L  complements   l[3];
#    g7: 0 <= y[4] <= L  complements   l[4];

using JuMP, Complementarity
using Ipopt
using Test


@testset "gnash1m.jl" begin

    gg = 5000^(1/g)

    gnash1m = Model(Ipopt.Optimizer)

    @variable(gnash1m, 0 <= x <= L)
    @variable(gnash1m, y[1:4])
    @variable(gnash1m, l[1:4])
    @variable(gnash1m, Q >= 0)
    @constraint(gnash1m, Q == x+y[1]+y[2]+y[3]+y[4])
    @NLobjective(gnash1m, Min, c[1]*x + b[1]/(b[1]+1)*K[1]^(-1/b[1])*x^((1+b[1])/b[1])
     		             - x*( gg*Q^(-1/g) ) )

    @NLconstraint(gnash1m, cnstr[i=1:4], 0 == ( c[i+1] + K[i+1]^(-1/b[i+1])*y[i] ) - ( gg*Q^(-1/g) )
                                      - y[i]*( -1/g*gg*Q^(-1-1/g) ) - l[i] )


    @macroexpand @complements(gnash1m, l[1], 0<=y[1]<=L,:smooth)
    
    
    
    for i in 1:4
        @complements(gnash1m, l[i], 0 <= y[i] <= L, :smooth)
    end

    # Initial solutions to help reaching the optimality
    set_start_value.(y, [11, 19, 20, 20])
    set_start_value.(l, [0, 0, -3, -5])
    set_start_value(x, 5)
    set_start_value(Q, 70)

    JuMP.optimize!(gnash1m)

    @show JuMP.value.(y)
    @show JuMP.value.(l)
    @show JuMP.value(x)
    @show JuMP.value(Q)

    # JuMP.value.(y) = [11.3327, 19.5003, 20.0, 20.0]
    # JuMP.value.(l) = [-1.35686e-10, -9.7489e-9, -2.51077, -5.18082]
    # JuMP.value(x) = 6.369341692612341
    # JuMP.value(Q) = 77.20235592852553
    # JuMP.objective_value(gnash1m) = -6.116708234438121

    @show JuMP.objective_value(gnash1m)
    @test isapprox(JuMP.objective_value(gnash1m), -6.11671, atol=1e-4)
    @test isapprox( JuMP.value.(l)' * (L .- JuMP.value.(y)), 0, atol=1e-5)
end
