type OneDimension
    lb::Float64
    var::JuMP.Variable
    ub::Float64
    F::JuMP.NonlinearExpression
end

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

function correspond(m::Model, F::JuMP.NonlinearExpression, var::JuMP.Variable)
    lb = getLower(var)
    ub = getUpper(var)
    new_dimension = OneDimension(float(lb), var, float(ub), F)
    data = getMCPData(m)
    push!(data, new_dimension)
end

# # These vectorinzed functions need more tests.
# # Seems that it should be JuMPArray in general.
# function correspond(m::Model, F::Array{JuMP.NonlinearExpression,1}, var::Array{JuMP.Variable,1})
#     @assert length(F) == length(var)
#     for i in 1:length(var)
#         correspond(m, F[i], var[i])
#     end
# end

function correspond(m::Model, F::Array{JuMP.NonlinearExpression}, var::Array{JuMP.Variable})
    vars = collect(var)
    Fs = collect(F)

    @assert length(vars) == length(Fs)

    for i in 1:length(vars)
        correspond(m, Fs[i], vars[i])
    end
end


# function correspond{T<:Any}(m::Model, F::JuMP.JuMPArray{JuMP.NonlinearExpression,1,Tuple{T}}, var::JuMP.JuMPArray{JuMP.Variable,1,Tuple{T}})
#
#     @assert length(var.innerArray) == length(F.innerArray)
#
#     for i in 1:length(var.innerArray)
#         correspond(m, F.innerArray[i], var.innerArray[i])
#     end
# end

function correspond(m::Model, F::JuMP.JuMPArray, var::JuMP.JuMPArray)
    vars = collect(var.innerArray)
    Fs = collect(F.innerArray)

    @assert length(vars) == length(Fs)

    for i in 1:length(vars)
        correspond(m, Fs[i], vars[i])
    end
end


# ::JuMP.JuMPArray{JuMP.NonlinearExpression,1,Tuple{UnitRange{Int64}}}, ::JuMP.JuMPArray{JuMP.Variable,1,Tuple{UnitRange{Int64}}})
#
# ::JuMP.JuMPArray{JuMP.NonlinearExpression,1,Tuple{Array{ASCIIString,1}}}, ::JuMP.JuMPArray{JuMP.Variable,1,Tuple{Array{ASCIIString,1}}})


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

        # This part might need to use MathProgBase function evaluations.
        # Just like below in myjac(z)
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


    # After solving set the values in m::JuMP.Model to the solution obtained.
    for i in 1:n
        variable = data[i].var
        value = z[getLinearIndex(variable)]
        setValue(variable, value)
    end

    # This function has changed the content of m already.
    return m
end
