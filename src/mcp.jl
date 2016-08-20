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



# The most basic one.
function complements(m::Model, F::JuMP.NonlinearExpression, var::JuMP.Variable)
    lb = getlowerbound(var)
    ub = getupperbound(var)
    new_dimension = ComplementarityType(lb, var, ub, F, linearindex(var))
    mcp_data = getMCPData(m)
    push!(mcp_data, new_dimension)
end
function complements(m::Model, var::JuMP.Variable, F::JuMP.NonlinearExpression)
    complements(m, F, var)
end



function complements(m::Model, F::Array{JuMP.NonlinearExpression}, var::Array{JuMP.Variable})
    vars = collect(var)
    Fs = collect(F)

    @assert length(vars) == length(Fs)

    for i in 1:length(vars)
        complements(m, Fs[i], vars[i])
    end
end
function complements(m::Model, var::Array{JuMP.Variable}, F::Array{JuMP.NonlinearExpression}, )
    complements(m::Model, F, var)
end

function complements(m::Model,  arr1::JuMP.JuMPArray, arr2::JuMP.JuMPArray)
    arr1s = collect(arr1.innerArray)
    arr2s = collect(arr2.innerArray)

    @assert length(arr1s) == length(arr2s)

    for i in 1:length(arr1s)
        complements(m, arr1s[i], arr2s[i])
    end
end




function getBoundsLinearIndex(mcp_data)
    n = length(mcp_data)
    lb = zeros(n)
    ub = ones(n)
    for i in 1:n
        lb[linearindex(mcp_data[i].var)] = mcp_data[i].lb
        ub[linearindex(mcp_data[i].var)] = mcp_data[i].ub
    end
    return lb, ub
end



# Placeholder for multiple methods in the future
function solveMCP(m::Model; method=:path)
    return _solve_path(m)
end

function sortMCPDataperm(obj::Array{ComplementarityType,1})
    n = length(obj)
    ref = Array{Int}(n)
    for i in 1:n
        ref[i] = obj[i].lin_idx
    end

    return sortperm(ref)
end

# Using PATHSolver
function _solve_path(m::Model)

    function myfunc(z)
        # z is in LindexIndex, passed from PATHSolver

        d = JuMP.NLPEvaluator(m)
        MathProgBase.initialize(d, [:Grad])
        F_val = zeros(n)
        MathProgBase.eval_g(d, F_val, z)

        # F_val also should be in LindexIndex
        # since it is the order in which constraints are added

        return F_val
    end

    function myjac(z)
        # z is in LindexIndex, passed from PATHSolver

        d = JuMP.NLPEvaluator(m)
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
    # Add constraint in the order of LinearIndex
    p = sortMCPDataperm(mcp_data)
    @NLconstraint(m, dummy[i=1:n], mcp_data[p[i]].F == 0)

    # lb and ub in LinearIndex
    lb, ub = getBoundsLinearIndex(mcp_data)

    # Solve the MCP using PATHSolver
    # ALL inputs to PATHSolver must be in LinearIndex
    status, z, f = PATHSolver.solveMCP(myfunc, myjac, lb, ub)
    # z, f are in LinearIndex

    # After solving set the values in m::JuMP.Model to the solution obtained.
    for i in 1:n
        setvalue(mcp_data[i].var, z[mcp_data[i].lin_idx])
    end

    # This function has changed the content of m already.
    return status
end




# The below will be useful for creating wrapper macros
# so that users don't need to do 'using Complementarity, JuMP'
# but, just 'using Complementarity'

# Modification of deprecate_macro from JuMP.jl
macro macro_wrapper(old, new)
    oldmac = symbol(string("@",old))
    newmac = symbol(string("@",new))
    s = string(oldmac," is just a wrapper of Jump.", newmac, ".")
    if VERSION > v"0.5-"
        # backtraces are ok on 0.5
        # depwarn = :(Base.depwarn($s,$(Meta.quot(oldmac))))
        depwarn = :()
    else
        # backtraces are junk on 0.4
        # depwarn = :(Base.warn_once($s))
        depwarn = :()
    end
    @eval macro $old(args...)
        return Expr(:block, $depwarn, Expr(:macrocall, $(Meta.quot(newmac)), [esc(x) for x in args]...))
    end
    eval(Expr(:export,oldmac))
    return
end

@macro_wrapper mapping defNLExpr
