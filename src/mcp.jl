mutable struct ComplementarityType
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
    m.ext[:MCP] = Array{ComplementarityType}(undef, 0)
    return m
end

function getMCPData(m::Model)
    if haskey(m.ext, :MCP)
        return m.ext[:MCP]::Array
    else
        error("The 'getMCPData' function is only for MCP models as in ComplementarityType.jl")
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
    var_name = Array{String}(undef, n)
    F_name = Array{String}(undef, n)
    for i in 1:n
        var_name[linearindex(mcp_data[i].var)] = mcp_data[i].var_name
        F_name[linearindex(mcp_data[i].var)] = mcp_data[i].F_name
    end
    return var_name, F_name
end

function getInitialValuesLinearIndex(m, mcp_data)
    n = length(mcp_data)
    initial_values = Array{Float64}(undef, n)
    for i in 1:n

        # init_val = getvalue(mcp_data[i].var) # This throws a warning when it is not set.
        init_val = m.colVal[linearindex(mcp_data[i].var)]
        if isnan(init_val)
            # When an initial value is not provided by 'setvalue(x, 1.0)'
            # it is set to the lower bound
            init_val = mcp_data[i].lb
        end
        initial_values[linearindex(mcp_data[i].var)] = init_val
    end
    return initial_values
end

# Placeholder for multiple methods in the future
function solveMCP(m::Model; solver=:PATH, method=:trust_region, linear=false)
    if solver == :PATH
        return _solve_path(m, linear=linear)
    elseif solver == :NLsolve
        return _solve_nlsolve(m, method=method)
    end
end

function solveLCP(m::Model; solver=:PATH, method=:trust_region)
    solveMCP(m, solver=solver, method=method, linear=true)
end

function sortMCPDataperm(obj::Array{ComplementarityType,1})
    n = length(obj)
    ref = Array{Int}(undef, n)
    for i in 1:n
        ref[i] = obj[i].lin_idx
    end

    return sortperm(ref)
end

# Using PATHSolver
function _solve_path(m::Model; linear=false)

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

    # Declaring MCP mapping F as constraints
    # in order to query Jacobian using AutoDiff thru MathProgBase
    # i = LinearIndex
    # Add constraint in the order of LinearIndex
    p = sortMCPDataperm(mcp_data)
    @NLconstraint(m, [i=1:n], mcp_data[p[i]].F == 0)

    # lb and ub in LinearIndex
    lb, ub = getBoundsLinearIndex(mcp_data)

    # var_name, F_name in LinearIndex
    var_name, F_name = getNamesLinearIndex(mcp_data)

    # initial values
    initial_values = getInitialValuesLinearIndex(m, mcp_data)

    # Solve the MCP using PATHSolver
    # ALL inputs to PATHSolver must be in LinearIndex
    if linear==true
        J0 = myjac(zeros(size(lb)))
        Jr = myjac(100*rand(Float64, size(ub)))

        if norm(J0-Jr, 1) > 10e-8
            error("The mappings do not seem linear. Rerun 'solveMCP()' after removing 'linear=true'.")
        end
        status, z, f = PATHSolver.solveLCP(myfunc, J0, lb, ub, initial_values, var_name, F_name)
    else
        status, z, f = PATHSolver.solveMCP(myfunc, myjac, lb, ub, initial_values, var_name, F_name)
    end

    # z, f are in LinearIndex

    # After solving set the values in m::JuMP.Model to the solution obtained.
    for i in 1:n
        setvalue(mcp_data[i].var, z[mcp_data[i].lin_idx])
    end



    # Cleanup. Remove all dummy @NLconstraints added,
    # so that the model can be re-used for multiple runs
    m.nlpdata.nlconstr = Array{JuMP.GenericRangeConstraint{JuMP.NonlinearExprData},1}(undef, 0)


    # This function has changed the content of m already.
    return status
end


function _solve_nlsolve(m::Model; method=:trust_region)

    function myfunc!(fvec, z)
        # z is in LindexIndex, passed from PATHSolver

        d = JuMP.NLPEvaluator(m)
        MathProgBase.initialize(d, [:Grad])
        F_val = zeros(n)
        MathProgBase.eval_g(d, F_val, z)

        copyto!(fvec, F_val)

        # F_val also should be in LindexIndex
        # since it is the order in which constraints are added
    end

    function myjac!(fjac, z)
        # z is in LindexIndex, passed from PATHSolver

        d = JuMP.NLPEvaluator(m)
        MathProgBase.initialize(d, [:Grad])
        I,J = MathProgBase.jac_structure(d)
        jac_val = zeros(size(J))
        MathProgBase.eval_jac_g(d, jac_val, z)

        # return matrix also should be in LindexIndex
        # since it is the order in which constraints are added

        sparse_fjac = sparse(I, J, jac_val)
        copy!(fjac, full(sparse_fjac))
    end

    mcp_data = getMCPData(m)
    n = length(mcp_data)

    # Two Indices
    # MCP_Index: the order stored in MCPModel = array index of Array{ComplementarityType}
    # LinearIndex: the order used in JuMP / MathProgBase

    # Declaring MCP mapping F as constraints
    # in order to query Jacobian using AutoDiff thru MathProgBase
    # i = LinearIndex
    # Add constraint in the order of LinearIndex
    p = sortMCPDataperm(mcp_data)
    @NLconstraint(m, [i=1:n], mcp_data[p[i]].F == 0)

    # lb and ub in LinearIndex
    lb, ub = getBoundsLinearIndex(mcp_data)

    # setting initial_z
    initial_z = (lb+ub) / 2
    for i = 1:length(initial_z)
        if initial_z[i] == Inf
            initial_z[i] = lb[i]
        end
    end

    # If the user has provided any initial_z, use it.
    for i in 1:n
        gv = JuMP._getValue(mcp_data[i].var)
        if !isnan(gv)
          initial_z[mcp_data[i].lin_idx] = gv
        end
    end

    # Solve the MCP using NLsolve
    # ALL inputs to NLsolve must be in LinearIndex

    r = NLsolve.mcpsolve(myfunc!, myjac!, lb, ub, initial_z, method = method,
        iterations = 10_000)
    # function mcpsolve{T}(f,
    #                   j,
    #                   lower::Vector,
    #                   upper::Vector,
    #                   initial_x::AbstractArray{T};
    #                   method::Symbol = :trust_region,
    #                   reformulation::Symbol = :smooth,
    #                   xtol::Real = zero(T),
    #                   ftol::Real = convert(T,1e-8),
    #                   iterations::Integer = 1_000,
    #                   store_trace::Bool = false,
    #                   show_trace::Bool = false,
    #                   extended_trace::Bool = false,
    #                   linesearch = LineSearches.BackTracking(),
    #                   factor::Real = one(T),
    #                   autoscale = true,
    #                   inplace = true)

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


    # Cleanup. Remove all dummy @NLconstraints added,
    # so that the model can be re-used for multiple runs
    m.nlpdata.nlconstr = Array{JuMP.GenericRangeConstraint{JuMP.NonlinearExprData},1}(undef, 0)



    # This function has changed the content of m already.
    return r
end






# https://stackoverflow.com/questions/50084877/how-to-alias-a-macro-in-julia#comment87277306_50085297
@eval const $(Symbol("@mapping")) = $(Symbol("@NLexpression"))

# macro mapping(args...)
#   if length(args) != 3
#     error("3 arguments are required in @mapping(...)")
#   end
#   quote
#     @NLexpression $(args...)
#   end
# end

# macro mapping(args...)
#   if length(args) != 3
#     error("3 arguments are required in @mapping(...)")
#   end
#
#   m = args[1]
#   name = args[2]
#   ex = args[3]
#
#   expression = Expr(:macrocall, Symbol("@NLexpression"), m, name, ex)
#
#   return esc(expression)
# end


function add_complementarity(m::JuMP.Model, var::JuMP.Variable, F::JuMP.NonlinearExpression, F_name::String)
  lb = getlowerbound(var)
  ub = getupperbound(var)
  var_name = getname(var)
  new_dimension = ComplementarityType(lb, var, ub, F, linearindex(var), var_name, F_name)
  mcp_data = getMCPData(m)
  push!(mcp_data, new_dimension)
end

macro complementarity(m, F, var)
  F_base_name = string(F)
  F_sym = string(F)
  var_sym = string(var)

  m = esc(m)
  F = esc(F)
  var = esc(var)

  quote

    if isa($F, JuMP.JuMPArray) || isa($F, Array)
      @assert length($F) == length($var)
    end

    # when var is a single JuMP variable
    if isa($var, JuMP.Variable)
      ex_var = Meta.parse(getname($var))
      ex_F = Meta.parse($F_base_name)

      if typeof(ex_var) == Symbol || typeof(ex_F) == Symbol
        add_complementarity($m, $var, $F, $F_base_name)

      elseif typeof(ex_var) == Expr && typeof(ex_F) == Expr
        @assert ex_var.head == :ref
        @assert ex_F.head == :ref

        ex_F.args[2] = ex_var.args[2]
        F_name = string(ex_F)

        add_complementarity($m, $var, $F, F_name)

      else
        error("Error in @complementarity. Please file an issue with an example.")
      end

    # when var is a single dimensional Array of JuMP.Variable
    elseif isa($var, Array{JuMP.Variable, 1})

      idx_list = 1:length($var)

      for idx in idx_list
        idx_name = idx
        F_name = string($F_base_name, "[", idx_name, "]")
        var_idx = $var[idx]
        F_idx = $F[idx]

        add_complementarity($m, var_idx, F_idx, F_name)
      end

    else # isa($var, JuMP.JuMPArray) && length(($var).indexsets) == 1 or > 1
      # when var is a multi-dimensional JuMP variable array
      ex = :(Base.product())
      for i in 1:length($var.indexsets)
        push!(ex.args, $var.indexsets[i])
      end
      idx_list = collect(eval(ex))

      for idx in idx_list
        F_name = string($F_base_name, "[", join(idx,","), "]")

        var_idx = $var[idx...]
        F_idx = $F[idx...]

        add_complementarity($m, var_idx, F_idx, F_name)
      end

    end # end of if


  end # end of quote



end # end of @complementarity





#
