type ComplementarityType
    lb::Float64
    var::JuMP.Variable
    ub::Float64
    F::JuMP.NonlinearExpression
    lin_idx::Int
end

function MCPModel()
    m = Model()
    m.ext[:MCP] = Array(ComplementarityType,0)
    return m
end

function getMCPData(m::Model)
    if haskey(m.ext, :MCP)
        return m.ext[:MCP]::Array
    else
        error("The 'getMCPData' function is only for MCP models as in ComplementarityType.jl")
    end
end

function correspond(m::Model, F::JuMP.NonlinearExpression, var::JuMP.Variable)
    lb = getLower(var)
    ub = getUpper(var)
    new_dimension = ComplementarityType(lb, var, ub, F, getLinearIndex(var))
    mcp_data = getMCPData(m)
    push!(mcp_data, new_dimension)
end

function correspond(m::Model, F::Array{JuMP.NonlinearExpression}, var::Array{JuMP.Variable})
    vars = collect(var)
    Fs = collect(F)

    @assert length(vars) == length(Fs)

    for i in 1:length(vars)
        correspond(m, Fs[i], vars[i])
    end
end

function correspond(m::Model, F::JuMP.JuMPArray, var::JuMP.JuMPArray)
    vars = collect(var.innerArray)
    Fs = collect(F.innerArray)

    @assert length(vars) == length(Fs)

    for i in 1:length(vars)
        correspond(m, Fs[i], vars[i])
    end
end

function getBoundsLinearIndex(mcp_data)
    n = length(mcp_data)
    lb = zeros(n)
    ub = ones(n)
    for i in 1:n
        lb[getLinearIndex(mcp_data[i].var)] = mcp_data[i].lb
        ub[getLinearIndex(mcp_data[i].var)] = mcp_data[i].ub
    end
    return lb, ub
end



# Placeholder for multiple methods in the future
function solveMCP(m::Model; method=:path)
    return _solve_path(m)
end



# Using PATHSolver
function _solve_path(m::Model)

    function myfunc0(z)
        mcp_data = getMCPData(m)
        F_ret = similar(z)

        for i in 1:length(mcp_data)
            setValue(mcp_data[i].var, z[i])
        end

        for i in 1:length(mcp_data)
            F_ret[i] = getValue(mcp_data[i].F)
        end

        return F_ret
    end

    function myfunc(z)
        # z is in LindexIndex, passed from PATHSolver

        d = JuMPNLPEvaluator(m)
        MathProgBase.initialize(d, [:Grad])
        F_val = zeros(n)
        MathProgBase.eval_g(d, F_val, z)

        # F_val also should be in LindexIndex
        # since it is the order in which constraints are added

        return F_val
    end

    function myjac(z)
        # z is in LindexIndex, passed from PATHSolver

        d = JuMPNLPEvaluator(m)
        MathProgBase.initialize(d, [:Grad])
        I,J = MathProgBase.jac_structure(d)
        jac_val = zeros(size(J))
        MathProgBase.eval_jac_g(d, jac_val, z)

        # return matrix also should be in LindexIndex
        # since it is the order in which constraints are added

        return sparse(I, J, jac_val)
    end

    mcp_data = getMCPData(m)
    n = length(mcp_data)

    # Two Indices
    # MCP_Index: the order stored in MCPModel = array index of Array{ComplementarityType}
    # LinearIndex: the order used in JuMP / MathProgBase

    # Declaring MCP operator F as constraints
    # in order to query Jacobian using AutoDiff thru MathProgBase
    # i = LinearIndex
    # Add constraint in the order of LinearIndex (NEEDS IMPROVEMENT)
    constraint = Array{Any}(n)
    for i in 1:n
        for j in 1:n
            if mcp_data[j].lin_idx == i
                constraint[i] = @addNLConstraint(m, mcp_data[j].F == 0)
            end
        end
    end
    # @addNLConstraint(m, constraint[i=1:n], mcp_data[i].F == 0)

    # lb and ub in LinearIndex
    lb, ub = getBoundsLinearIndex(mcp_data)


    # Solve the MCP using PATHSolver
    # ALL inputs to PATHSolver must be in LinearIndex
    z, f = PATHSolver.solveMCP(myfunc, myjac, lb, ub)
    # z, f are in LinearIndex

    # After solving set the values in m::JuMP.Model to the solution obtained.
    for i in 1:n
        setValue(mcp_data[i].var, z[mcp_data[i].lin_idx])
    end

    # This function has changed the content of m already.
    return m
end
