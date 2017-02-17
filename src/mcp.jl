type ComplementarityType
    lb::Float64
    var::JuMP.Variable
    ub::Float64
    F::JuMP.NonlinearExpression
    lin_idx::Int
    var_name::String
    F_name::String
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
    Base.depwarn("complements is deprecated. Use @complementarity instead.", :complements)
    lb = getlowerbound(var)
    ub = getupperbound(var)
    var_name = getname(var)
    F_name = ""
    new_dimension = ComplementarityType(lb, var, ub, F, linearindex(var), var_name, F_name)
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


function getNamesLinearIndex(mcp_data)
    n = length(mcp_data)
    var_name = Array{String}(n)
    F_name = Array{String}(n)
    for i in 1:n
        var_name[linearindex(mcp_data[i].var)] = mcp_data[i].var_name
        F_name[linearindex(mcp_data[i].var)] = mcp_data[i].F_name
    end
    return var_name, F_name
end


# Placeholder for multiple methods in the future
function solveMCP(m::Model; solver=:PATH, method=:trust_region)
    if solver == :PATH
        return _solve_path(m)
    elseif solver == :NLsolve
        return _solve_nlsolve(m, method=method)
    end
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

    # var_name, F_name in LinearIndex
    var_name, F_name = getNamesLinearIndex(mcp_data)

    # Solve the MCP using PATHSolver
    # ALL inputs to PATHSolver must be in LinearIndex
    status, z, f = PATHSolver.solveMCP(myfunc, myjac, lb, ub, var_name, F_name)
    # z, f are in LinearIndex

    # After solving set the values in m::JuMP.Model to the solution obtained.
    for i in 1:n
        setvalue(mcp_data[i].var, z[mcp_data[i].lin_idx])
    end

    # This function has changed the content of m already.
    return status
end


function _solve_nlsolve(m::Model; method=:trust_region)

    function myfunc!(z, fvec)
        # z is in LindexIndex, passed from PATHSolver

        d = JuMP.NLPEvaluator(m)
        MathProgBase.initialize(d, [:Grad])
        # F_val = zeros(n)
        MathProgBase.eval_g(d, fvec, z)

        # F_val also should be in LindexIndex
        # since it is the order in which constraints are added

        # return fvec
    end

    function myjac!(z, fjac)
        # z is in LindexIndex, passed from PATHSolver

        d = JuMP.NLPEvaluator(m)
        MathProgBase.initialize(d, [:Grad])
        I,J = MathProgBase.jac_structure(d)
        jac_val = zeros(size(J))
        MathProgBase.eval_jac_g(d, jac_val, z)

        # return matrix also should be in LindexIndex
        # since it is the order in which constraints are added

        sparse_fjac = sparse(I, J, jac_val)
        fjac = full(sparse_fjac)
        # return fjac
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
    initial_z = (lb+ub) / 2
    for i = 1:length(initial_z)
        if initial_z[i] == Inf
            initial_z[i] = lb[i]
        end
    end

    # Solve the MCP using PATHSolver
    # ALL inputs to PATHSolver must be in LinearIndex
    # status, z, f = PATHSolver.solveMCP(myfunc, myjac, lb, ub)

    r = NLsolve.mcpsolve(myfunc!, myjac!, lb, ub, initial_z, method = method,
        iterations = 10_000)
    # function mcpsolve{T}(f!::Function,
    #                   g!::Function,
    #                   lower::Vector,
    #                   upper::Vector,
    #                   initial_x::Vector{T};
    #                   method::Symbol = :trust_region,
    #                   reformulation::Symbol = :smooth,
    #                   xtol::Real = zero(T),
    #                   ftol::Real = convert(T,1e-8),
    #                   iterations::Integer = 1_000,
    #                   store_trace::Bool = false,
    #                   show_trace::Bool = false,
    #                   extended_trace::Bool = false,
    #                   linesearch!::Function = Optim.backtracking_linesearch!,
    #                   factor::Real = one(T),
    #                   autoscale::Bool = true)


    # julia> fieldnames(r)
    # 12-element Array{Symbol,1}:
    #  :method
    #  :initial_x
    #  :zero
    #  :residual_norm
    #  :iterations
    #  :x_converged
    #  :xtol
    #  :f_converged
    #  :ftol
    #  :trace
    #  :f_calls
    #  :g_calls


    # r.zero is in LinearIndex

    # After solving set the values in m::JuMP.Model to the solution obtained.
    for i in 1:n
        setvalue(mcp_data[i].var, r.zero[mcp_data[i].lin_idx])
    end

    # This function has changed the content of m already.
    return r
end









macro operator(args...)
  if length(args) != 3
    error("3 arguments are required in @operator(...)")
  end

  m = args[1]
  name = args[2]
  ex = args[3]

  expression = Expr(:macrocall, Symbol("@NLexpression"), m, name, ex)

  return esc(expression)
end



macro complementarity(m, F, var)
  F_base_name = string(F)
  F_sym = string(F)
  var_sym = string(var)

  m = esc(m)
  F = esc(F)
  var = esc(var)

  quote
    if isa($F, JuMP.JuMPArray)
      @assert length($F) == length($var)
    end

    if isa($var, JuMP.Variable)
      # when var is a single JuMP variable
      lb = getlowerbound($var)
      ub = getupperbound($var)
      var_name = getname($var)
      F_name = $F_base_name
      new_dimension = ComplementarityType(lb, $var, ub, $F, linearindex($var), var_name, F_name)
      mcp_data = getMCPData($m)
      push!(mcp_data, new_dimension)

    elseif isa($var, Array{JuMP.Variable, 1})
      # when var is a single dimensional Array of JuMP.Variable
      idx_list = 1:length($var)

      for idx in idx_list
        idx_name = idx
        F_name = string($F_base_name, "[", idx_name, "]")
        var_idx = $var[idx]
        F_idx = $F[idx]

        lb = getlowerbound(var_idx)
        ub = getupperbound(var_idx)
        var_name = getname(var_idx)
        new_dimension = ComplementarityType(lb, var_idx, ub, F_idx, linearindex(var_idx), var_name, F_name)
        mcp_data = getMCPData($m)
        push!(mcp_data, new_dimension)
      end

    elseif isa($var, JuMP.JuMPArray) && length(($var).indexsets) == 1
      # when var is a single dimensional JuMPArray of JuMP.Variable
      idx_list = $var.indexsets[1]
      for idx in idx_list
        idx_name = idx
        F_name = string($F_base_name, "[", idx_name, "]")
        var_idx = $var[idx]
        F_idx = $F[idx]

        lb = getlowerbound(var_idx)
        ub = getupperbound(var_idx)
        var_name = getname(var_idx)
        new_dimension = ComplementarityType(lb, var_idx, ub, F_idx, linearindex(var_idx), var_name, F_name)
        mcp_data = getMCPData($m)
        push!(mcp_data, new_dimension)
      end

    else
      # when var is a multi-dimensional JuMP variable array
      ex = :(Base.product())
      for i in 1:length($var.indexsets)
        push!(ex.args, $var.indexsets[i])
      end
      idx_list = collect(eval(ex))

      for idx in idx_list
        idx_name = idx
        F_name = string($F_base_name, "[", idx_name, "]")
        F_name = replace(F_name, "\"", "")
        F_name = replace(F_name, "(", "")
        F_name = replace(F_name, ")", "")

        F_name = string($F_base_name, "[", join(idx_name,","), "]")

        var_idx = $var[idx...]
        F_idx = $F[idx...]

        # if length(idx) == 1
        #   var_idx = $var[idx]
        #   F_idx = $F[idx]
        # else
        #   var_ex = Expr(:ref, Symbol($var_sym))
        #   F_ex = Expr(:ref, Symbol($F_sym))
        #   for j in 1:length(idx)
        #     push!(var_ex.args, idx[j])
        #     push!(F_ex.args, idx[j])
        #   end
        #   var_idx = var_ex
        #   F_idx = F_ex
        # end
        lb = getlowerbound(var_idx)
        ub = getupperbound(var_idx)
        var_name = getname(var_idx)
        new_dimension = ComplementarityType(lb, var_idx, ub, F_idx, linearindex(var_idx), var_name, F_name)
        mcp_data = getMCPData($m)
        push!(mcp_data, new_dimension)
      end
    end
  end # end of quote



end # end of @complementarity





#
