
esc_nonconstant(x::Number) = x
esc_nonconstant(x) = esc(x)
complements_error(args, str) = error("In @complements($(join(args,","))): ", str)

# for backward compatibility with 0.4
function comparison_to_call(ex)
    if isexpr(ex,:comparison) && length(ex.args) == 3
        return Expr(:call, ex.args[2], ex.args[1], ex.args[3])
    else
        return ex
    end
end


function smooth(c1, c2)
    c1s = Expr(:call, esc(:(^)), c1, 2)
    c2s = Expr(:call, esc(:(^)), c2, 2)
    cs = Expr(:call, esc(:(+)), c1s, c2s, :($mpec_tol))
    csqrt = Expr(:call, esc(:sqrt), cs)

    csum = Expr(:call, esc(:+), c1, c2)
    cc = Expr(:call, esc(:-), csqrt, csum)

    return Expr(:call, esc(:(==)), cc, 0  )
end

function get_complementarity(c1, c2, method)

    if method == :smooth
        expr = smooth(c1, c2)
    elseif method == :simple
        cc = Expr(:call, esc(:(*)), c1, c2)
        expr = Expr(:call, esc(:(<=)), cc, :($mpec_tol)  )
    else
        expr = smooth(c1, c2)
    end

    return expr

end




macro complements(args...)

    # Currently supports
    #  F ⟂ lb <= x <= ub
    #  or
    #  0 <= F ⟂ x >= 0
    #  0 >= F ⟂ x >= 0
    #  0 <= F ⟂ x <= 0
    #  0 >= F ⟂ x <= 0
    #  Note: the bound on F must be on the left hand side, and
    #        the obund on x must be on the right hand side.

    # Many parts of this macro were copied from @variable of JuMP.jl
    # Thanks!

    m = esc(args[1])
    F = args[2]
    x = args[3]

    method = :smooth
    mpec_tol = :(1e-8)
    if length(args) > 3
        method = args[4]
    end

    ub_F = Inf
    lb_F = -Inf
    ub_x = Inf
    lb_x = -Inf

    xhaslb = false
    xhasub = false
    Fhaslb = false
    Fhasub = false


    # Identify the variable bounds. Five (legal) possibilities are "x >= lb",
    # "x <= ub", "lb <= x <= ub", "x == val", or just plain "x"
    if VERSION < v"0.5.0-dev+3231"
        x = comparison_to_call(x)
        F = comparison_to_call(F)
    end


    ############################### x      #####################################
    if isexpr(x, :comparison) # two-sided

        if x.args[2] == :>= || x.args[2] == :≥
            # ub >= x >= lb
            x.args[4] == :>= || x.args[4] == :≥ || complements_error(args, "Invalid variable bounds")
            var = x.args[3]
            lb_x = esc_nonconstant(x.args[5])
            ub_x = esc_nonconstant(x.args[1])

            xhaslb = true
            xhasub = true
        elseif x.args[2] == :<= || x.args[2] == :≤
            # lb <= x <= ub
            var = x.args[3]
            (x.args[4] != :<= && x.args[4] != :≤) &&
                complements_error(args, "Expected <= operator after variable name.")

            lb_x = esc_nonconstant(x.args[1])
            ub_x = esc_nonconstant(x.args[5])

            xhaslb = true
            xhasub = true
        else
            complements_error(args, "Use the form lb <= ... <= ub.")
        end

    elseif isexpr(x, :call)
        if x.args[1] == :>= || x.args[1] == :≥
            # x >= lb
            var = x.args[2]
            @assert length(x.args) == 3
            lb_x = esc_nonconstant(x.args[3])
            ub_x = Inf
            xhaslb = true
            xhasub = false

            # # May also be ub >= x
            # if isa(eval(var), Number)
            #     var = x.args[3]
            #     @assert length(x.args) == 3
            #     ub_x = esc_nonconstant(x.args[2])
            #     lb_x = -Inf
            #     xhaslb = false
            #     xhasub = true
            # end

        elseif x.args[1] == :<= || x.args[1] == :≤
            # x <= ub
            var = x.args[2]
            @assert length(x.args) == 3
            ub_x = esc_nonconstant(x.args[3])
            lb_x = -Inf
            xhaslb = false
            xhasub = true

            # # May also be lb <= x
            # if isa(eval(var), Number)
            #     var = x.args[3]
            #     @assert length(x.args) == 3
            #     ub_x = Inf
            #     lb_x = esc_nonconstant(x.args[2])
            #     xhaslb = true
            #     xhasub = false
            # end

        else
            # # Its a comparsion, but not using <= ... <=
            # complements_error(args, "Unexpected syntax $(string(x)).")
            # No bounds provided - free variable
            # If it isn't, e.g. something odd like f(x), we'll handle later
            var = x
            lb_x = -Inf
            ub_x = Inf
        end

    else
        # No bounds provided - free variable
        # If it isn't, e.g. something odd like f(x), we'll handle later
        var = x
        lb_x = -Inf
        ub_x = Inf
    end
    ############################### x      #####################################




    ############################### F      #####################################
    if isexpr(F, :comparison) # two-sided

        if F.args[2] == :>= || F.args[2] == :≥
            # ub >= x >= lb
            F.args[4] == :>= || F.args[4] == :≥ || complements_error(args, "Invalid function bounds")
            func = F.args[3]

            lb_F = esc_nonconstant(F.args[5])
            ub_F = esc_nonconstant(F.args[1])

            Fhaslb = true
            Fhasub = true
        elseif F.args[2] == :<= || F.args[2] == :≤
            # lb <= x <= ub
            func = F.args[3]
            (F.args[4] != :<= && F.args[4] != :≤) &&
                complements_error(args, "Expected <= operator after function expression.")
            lb_F = esc_nonconstant(F.args[1])
            ub_F = esc_nonconstant(F.args[5])

            Fhaslb = true
            Fhasub = true
        else
            complements_error(args, "Use the form lb <= ... <= ub.")
        end

    elseif isexpr(F, :call)
        if F.args[1] == :>= || F.args[1] == :≥
            # x >= lb
            # func = F.args[2]
            # @assert length(F.args) == 3
            # lb_F = esc_nonconstant(F.args[3])
            # ub_F = Inf
            # Fhaslb = true

            # May also be ub >= x
            # if isa(eval(func), Number)
            func = F.args[3]
            @assert length(F.args) == 3
            ub_F = esc_nonconstant(F.args[2])
            lb_F = -Inf
            Fhaslb = false
            Fhasub = true
            # end

        elseif F.args[1] == :<= || F.args[1] == :≤
            # # x <= ub
            # func = F.args[2]
            # @assert length(F.args) == 3
            # ub_F = esc_nonconstant(F.args[3])
            # lb_F = -Inf
            # Fhasub = true

            # May also be lb <= x
            # if isa(eval(func), Number)
            func = F.args[3]
            @assert length(F.args) == 3
            ub_F = Inf
            lb_F = esc_nonconstant(F.args[2])
            Fhaslb = true
            Fhasub = false
            # end

        else
            # Its a comparsion, but not using <= ... <=
            # complements_error(args, "Unexpected syntax $(string(F)).")
            # No bounds provided - free funciable
            # If it isn't, e.g. something odd like f(x), we'll handle later
            func = F
            lb_F = -Inf
            ub_F = Inf
        end

    else
        # No bounds provided - free funciable
        # If it isn't, e.g. something odd like f(x), we'll handle later
        func = F
        lb_F = -Inf
        ub_F = Inf
    end
    ############################### F      #####################################

    number_bounds = xhaslb + xhasub + Fhaslb + Fhasub
    if number_bounds != 2
        complements_error(args, "The total number of bounds on the function and the variable must be exactly two.")
    end

    var_esc = esc(var)
    func_esc = esc(func)

    if isexpr(var, :call)
        complements_error(args, "The second argument must be a variable, not an expression.")
    end

    code = quote end

    # Bound constraints are NOT added in :smooth method.
    if method == :simple || method == :smooth
        if Fhaslb
            expr = Expr(:call, esc(:(>=)), esc(func), :($lb_F)  )
            push!(code.args, Expr(:macrocall, Symbol("@NLconstraint"), m, expr ) )
        end
        if Fhasub
            expr = Expr(:call, esc(:(<=)), esc(func), :($ub_F)  )
            push!(code.args, Expr(:macrocall, Symbol("@NLconstraint"), m, expr ) )
        end
        if xhaslb
            expr = Expr(:call, esc(:(>=)), esc(var), :($lb_x)  )
            push!(code.args, Expr(:macrocall, Symbol("@NLconstraint"), m, expr ) )
        end
        if xhasub
            expr = Expr(:call, esc(:(<=)), esc(var), :($ub_x)  )
            push!(code.args, Expr(:macrocall, Symbol("@NLconstraint"), m, expr ) )
        end
    end

    # There must be a better way of writing the codes below...
    if Fhaslb && Fhasub
        complements_error(args, "Both bounds on the function expression is currently not supported.")

    elseif xhaslb && xhasub

        # Additional variables are defined
        # v defined
        v = copy(var)
        v.args[1] = Symbol(join(["v", gensym(), var.args[1]]))
        v_expr = Expr(:call, esc(:>=), esc(v), 0)
        push!(code.args, Expr(:macrocall, Symbol("@variable"), m, v_expr))

        # w defined
        w = copy(var)
        w.args[1] = Symbol(join(["w", gensym(), var.args[1]]))
        w_expr = Expr(:call, esc(:>=), esc(w), 0)
        push!(code.args, Expr(:macrocall, Symbol("@variable"), m, w_expr))

        # v - w = func
        vwdiff = Expr(:call, esc(:(-)), esc(v), esc(w))
        # expr = Expr(:comparison, vwdiff, esc(:(==)), esc(func))
        expr = Expr(:call, esc(:(==)), vwdiff, esc(func))
        push!(code.args, Expr(:macrocall, Symbol("@constraint"), m, expr))

        # v * (x-lb) = 0
        c1 = Expr(:call, esc(:(-)), esc(var), :($lb_x))
        expr = get_complementarity(c1, esc(v), method)
        push!(code.args, Expr(:macrocall, Symbol("@NLconstraint"), m, expr))

        # w * (ub - x) = 0
        c1 = Expr(:call, esc(:(-)), :($ub_x), esc(var))
        expr = get_complementarity(c1, esc(w), method)
        push!(code.args, Expr(:macrocall, Symbol("@NLconstraint"), m, expr))

    elseif Fhaslb && xhaslb
        c1 = Expr(:call, esc(:(-)), esc(var), :($lb_x))
        c2 = Expr(:call, esc(:(-)), esc(func), :($lb_F))
        expr = get_complementarity(c1, c2, method)
        push!(code.args, Expr(:macrocall, Symbol("@NLconstraint"), m, expr ) )

    elseif Fhaslb && xhasub
        c1 = Expr(:call, esc(:(-)), :($ub_x), esc(var))
        c2 = Expr(:call, esc(:(-)), esc(func), :($lb_F))
        expr = get_complementarity(c1, c2, method)
        push!(code.args, Expr(:macrocall, Symbol("@NLconstraint"), m, expr ) )

    elseif Fhasub && xhaslb
        c1 = Expr(:call, esc(:(-)), esc(var), :($lb_x))
        c2 = Expr(:call, esc(:(-)), :($ub_F), esc(func))
        expr = get_complementarity(c1, c2, method)
        push!(code.args, Expr(:macrocall, Symbol("@NLconstraint"), m, expr ) )

    elseif Fhasub && xhasub
        c1 = Expr(:call, esc(:(-)), :($ub_x), esc(var))
        c2 = Expr(:call, esc(:(-)), :($ub_F), esc(func))
        expr = get_complementarity(c1, c2, method)
        push!(code.args, Expr(:macrocall, Symbol("@NLconstraint"), m, expr ) )

    else
        complements_error(args, "NO VALID CASE")
    end

    return code

end
