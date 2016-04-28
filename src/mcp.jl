type OneDimension
    lb::Float64
    var::JuMP.Variable
    ub::Float64
    F::JuMP.NonlinearExpression
end

# type MCPData
#     dimensions::Array{OneDimension}
# end

function MCPModel()
    m = Model()
    m.ext[:MCP] = Array(OneDimension,0)
    return m
end

function getMCPData(m::Model)
    if haskey(m.ext, :MCP)
        return m.ext[:MCP]::Array
    else
        error("The 'getMCPData' function is only for MCP models as in Complementarity.jl")
    end
end


function correspond(m::Model, lb::Number, var::JuMP.Variable, ub::Number, F::JuMP.NonlinearExpression)
    new_dimension = OneDimension(float(lb), var, float(ub), F)
    data = getMCPData(m)
    push!(data, new_dimension)
end

# These vectorinzed functions need more tests.
# Seems that it should be JuMPArray in general.
function correspond{T<:Number}(m::Model, lb::Array{T, 1}, var::Array{JuMP.Variable,1}, ub::Array{T, 1}, F::Array{JuMP.NonlinearExpression,1})
    @assert length(F) == length(var) == length(lb) == length(ub)
    for i in 1:length(var)
        correspond(m, lb[i], var[i], ub[i], F[i])
    end
end

function correspond{T<:Number}(m::Model, lb::Array{T}, var::Array{JuMP.Variable}, ub::Array{T}, F::Array{JuMP.NonlinearExpression})
    vars = collect(var)
    Fs = collect(F)
    lbs = collect(lb)
    ubs = collect(ub)

    @assert length(vars) == length(Fs) == length(lbs) == length(ubs)

    for i in 1:length(vars)
        correspond(m, lbs[i], vars[i], ubs[i], Fs[i])
    end
end



function solveMCP(m::Model; method=:path)
    return _solve_path(m)
end



function _solve_path(m::Model)

    function myfunc(z)
        data = getMCPData(m)
        # println("length of z = ", length(z))
        # println("length(data) = ", length(data))

        F_ret = similar(z)

        for i in 1:length(data)
            # println(data[i].var)
            # println(z[i])
            setValue(data[i].var, z[i])
        end

        for i in 1:length(data)
            F_ret[i] = getValue(data[i].F)
        end

        return F_ret
    end

    function myjac(z)
        n = length(z)
        z_val = zeros(n)
        for i in 1:n
            z_val[getLinearIndex(data[i].var)] = z[i]
        end

        d = JuMPNLPEvaluator(m)
        MathProgBase.initialize(d, [:Grad])
        I,J = MathProgBase.jac_structure(d)

        jac_val = zeros(size(J))
        MathProgBase.eval_jac_g(d, jac_val, z_val)


        return sparse(I, J, jac_val)

        # jac_mat = zeros(n,n)
        # for i in 1:length(I)
        #     jac_mat[I[i], J[i]] = jac_val[i]
        # end
        # return sparse(jac_mat)
    end

    data = getMCPData(m)
    n = length(data)

    # Declaring MCP operator F as constraints
    # in order to query Jacobian using AutoDiff thru MathProgBase
    @addNLConstraint(m, constraint[i=1:n], data[i].F == 0)

    lb = zeros(n)
    ub = zeros(n)

    for i in 1:n
        lb[i] = data[i].lb
        ub[i] = data[i].ub
    end

    # Solve the MCP using PATHSolver
    z, f = PATHSolver.solveMCP(myfunc, myjac, lb, ub)

    for i in 1:n
        variable = data[i].var
        value = z[getLinearIndex(variable)]
        setValue(variable, value)
    end

    return m

end
