using JuMP
using Ipopt

macro my_test(m,F)
    return :(@NLconstraint($(m), $(esc(F)) >= 0))
end


m = Model(Ipopt.Optimizer)
@variable(m, x)
@my_test(m,x)
