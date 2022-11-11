# To escape variable only
# (x+y)^2 ==> (esc(:x)+esc(:y))^2, instead of esc(:((x+y)^2))
# y[i] ==> esc(:y)[esc(:i)]
# In generator like below, i is placed in parameters and i is not escaped.
# sum(z[i] for i in 1:10) ==> sum(esc(:z)[i] for i in 1:10)


esc_variable(ex, parameters=Symbol[]) = ex

function esc_variable(ex::Symbol, parameters=Symbol[])
    if ex in parameters
        return ex
    else
        return esc(ex)
    end
end

function esc_variable(ex::Expr, parameters=Symbol[])

    ex2 = copy(ex)
    if ex2.head == :call
        for i in 2:length(ex2.args)
            ex2.args[i] = esc_variable(ex2.args[i], parameters)
        end
    elseif ex2.head == :generator
        ex2.args[1] = esc_variable(ex2.args[1], parameters)
        for i in 2:length(ex2.args)
            push!(parameters, ex2.args[i].args[1])
        end
    elseif ex2.head == :ref
        for i in 1:length(ex2.args)
            ex2.args[i] = esc_variable(ex2.args[i], parameters)
        end
    elseif ex2.head == :escape
        # do nothing
    else
        @show ex2
        dump(ex2)
        @error("In esc_variable(ex): ex2.head == $(ex2.head). Error. Not supported. Report to: https://github.com/chkwon/Complementarity.jl/issues")
    end
    return ex2
end



esc_nonconstant(x::Number) = x
esc_nonconstant(x) = esc(x)
complements_error(args, str) = error("In @complements($(join(args,","))): ", str)



function smooth(c1, c2)
    c1s = Expr(:call, :(^), c1, 2)
    c2s = Expr(:call, :(^), c2, 2)
    cs = Expr(:call, :(+), c1s, c2s, :($mpec_tol))
    csqrt = Expr(:call, :sqrt, cs)

    csum = Expr(:call, :+, c1, c2)
    cc = Expr(:call, :-, csqrt, csum)

    return Expr(:call, :(==), cc, 0  )
end

function get_complementarity(c1, c2, method)
    c1 = esc_variable(c1)
    c2 = esc_variable(c2)
    if method == :smooth
        expr = smooth(c1, c2)
    elseif method == :simple
        cc = Expr(:call, :(*), c1, c2)
        expr = Expr(:call, :(<=), cc, :($mpec_tol)  )
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
                complements_error(args, "Expected <= mapping after variable name.")

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

        elseif x.args[1] == :<= || x.args[1] == :≤
            # x <= ub
            var = x.args[2]
            @assert length(x.args) == 3
            ub_x = esc_nonconstant(x.args[3])
            lb_x = -Inf
            xhaslb = false
            xhasub = true

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
                complements_error(args, "Expected <= mapping after function expression.")
            lb_F = esc_nonconstant(F.args[1])
            ub_F = esc_nonconstant(F.args[5])

            Fhaslb = true
            Fhasub = true
        else
            complements_error(args, "Use the form lb <= ... <= ub.")
        end

    elseif isexpr(F, :call)
        if F.args[1] == :>= || F.args[1] == :≥
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

    if isexpr(var, :call)
        complements_error(args, "The second argument must be a variable, not an expression.")
    end

    code = quote end

    # Bound constraints are NOT added in :smooth method.
    if method == :simple || method == :smooth

        if Fhaslb
            push!(code.args, quote
                $add_nonlinear_constraint($m, :($$(esc_variable(func)) >= $$(lb_F)))
                #@NLconstraint( $(m), $(esc_variable(func)) >= $(lb_F) )
            end )
        end
        if Fhasub
            push!(code.args, quote
                $add_nonlinear_constraint($m, :($$(esc_variable(func)) <= $$(ub_F)))
                #@NLconstraint( $(m), $(esc_variable(func)) <= $(ub_F) )
            end )
        end
        if xhaslb
            push!(code.args, quote
                $add_nonlinear_constraint($m, :($$(esc_variable(var)) >= $$(lb_x)))
                #@NLconstraint( $(m), $(esc_variable(var)) >= $(lb_x) )
            end )
        end
        if xhasub
            push!(code.args, quote
                $add_nonlinear_constraint($m, :($$(esc_variable(var)) >= $$(ub_x)))
                #@NLconstraint( $(m), $(esc_variable(var)) <= $(ub_x) )
            end )
        end
    end



    var = esc_variable(var)
    func = esc_variable(func)

    # # There must be a better way of writing the codes below...
    if Fhaslb && Fhasub
        complements_error(args, "Both bounds on the function expression is currently not supported.")

    elseif xhaslb && xhasub
        # v defined
        push!(code.args, quote
            $(esc(:v)) = @variable($(m), lower_bound=0)
        end)

        # w defined
        push!(code.args, quote
            $(esc(:w)) = @variable($(m), lower_bound=0)
        end)

        # v - w = func
        push!(code.args, quote
            $add_nonlinear_constraint($m, :( $$(esc(:v)) - $$(esc(:w)) ==  $$(esc_variable(func))))
            #@NLconstraint( $(m), $(esc(:v))-$(esc(:w)) == $(esc_variable(func)) )
        end )

        # v * (x - lb) = 0
        c1 = Expr(:call, :(-), var, lb_x)
        expr = :($get_complementarity($c1, $(esc(:v)), $method))
        push!(code.args, quote
            $add_nonlinear_constraint($m, $expr)
        end )

        # w * (ub - x) = 0
        c1 = Expr(:call, :(-), ub_x, var)
        expr = :($get_complementarity($c1, $(esc(:w)), $method))
        push!(code.args, quote
            $add_nonlinear_constraint($m, $expr)
        end )

    elseif Fhaslb && xhaslb
        c1 = Expr(:call, :(-), var, lb_x)
        c2 = Expr(:call, :(-), func, lb_F)
        expr = :($get_complementarity($c1,$c2, $method))
        push!(code.args, quote
            $add_nonlinear_constraint($m, $expr)
        end )

    elseif Fhaslb && xhasub
        c1 = Expr(:call, :(-), ub_x, var)
        c2 = Expr(:call, :(-), func, lb_F)
        expr = :($get_complementarity($c1, $c2, $method))
        push!(code.args, quote
            $add_nonlinear_constraint($m, $expr)
        end )

    elseif Fhasub && xhaslb
        c1 = Expr(:call, :(-), var, lb_x)
        c2 = Expr(:call, :(-), ub_F, func)
        expr = :($get_complementarity($c1, $c2, $method))
        push!(code.args, quote
            $add_nonlinear_constraint($m, $expr)
        end )
    elseif Fhasub && xhasub
        c1 = Expr(:call, :(-), ub_x, var)
        c2 = Expr(:call, :(-), ub_F, func)
        expr = :($get_complementarity($c1, $c2, $method))
        push!(code.args, quote
            $add_nonlinear_constraint($m, $expr)
        end )

    else
        complements_error(args, "NO VALID CASE")
    end



    return code
end
