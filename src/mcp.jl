mutable struct ComplementarityType
    lb::AbstractFloat
    var::JuMP.VariableRef
    ub::AbstractFloat
    F::JuMP.NonlinearExpression
    raw_idx::Integer
    var_name::String
    F_name::String
    result_value::AbstractFloat
end

function MCPModel()
    m = JuMP.Model()
    m.ext[:MCP] = ComplementarityType[]
    return m
end

function get_MCP_data(m::JuMP.Model)
    if haskey(m.ext, :MCP)
        return m.ext[:MCP]::Array
    else
        error("The 'get_MCP_data' function is only for MCP models as in ComplementarityType.jl")
    end
end

raw_index(v::JuMP.VariableRef) = JuMP.index(v).value

function get_bounds_in_raw_index(mcp_data)
    n = length(mcp_data)
    lb = zeros(n)
    ub = ones(n)
    for i in 1:n
        lb[raw_index(mcp_data[i].var)] = mcp_data[i].lb
        ub[raw_index(mcp_data[i].var)] = mcp_data[i].ub
    end
    return lb, ub
end


function get_names_in_raw_index(mcp_data)
    n = length(mcp_data)
    var_name = Array{String}(undef, n)
    F_name = Array{String}(undef, n)
    for i in 1:n
        var_name[raw_index(mcp_data[i].var)] = mcp_data[i].var_name
        F_name[raw_index(mcp_data[i].var)] = mcp_data[i].F_name
    end
    return var_name, F_name
end

function get_initial_values_in_raw_index(m, mcp_data)
    n = length(mcp_data)
    initial_values = Array{Float64}(undef, n)
    for i in 1:n
        init_val = JuMP.start_value(mcp_data[i].var)
        if init_val == nothing
            init_val = mcp_data[i].lb
        end
        initial_values[raw_index(mcp_data[i].var)] = init_val
    end
    return initial_values
end

function solveMCP!(m::JuMP.Model; solver=:PATH, method=:trust_region, kwargs...)
    if solver == :PATH
        return _solve_path!(m; kwargs...)
    elseif solver == :NLsolve
        return _solve_nlsolve!(m, method=method)
    end
end


function sortperm_MCP_data(obj::Array{ComplementarityType,1})
    n = length(obj)
    ref = Array{Int}(undef, n)
    for i in 1:n
        ref[i] = obj[i].raw_idx
    end

    return sortperm(ref)
end

# Using PATHSolver
function _solve_path!(m::JuMP.Model; kwargs...)
    d = JuMP.NLPEvaluator(m)

    function function_callback(n::Cint, z::Vector{Cdouble}, F_val::Vector{Cdouble})
        @assert n == length(z) == length(F_val)
        # z is in RawIndex, passed from PATHSolver
        MOI.initialize(d, [:Jac])
        MOI.eval_constraint(d, F_val, z)
        # F_val also should be in RawIndex
        # since it is the order in which constraints are added
        return Cint(0)
    end

    function j_eval(z::Vector{Cdouble})
        MOI.initialize(d, [:Jac])
        J_struct = MOI.jacobian_structure(d)
        I = first.(J_struct)
        J = last.(J_struct)

        jac_val = zeros(length(J))
        MOI.eval_constraint_jacobian(d, jac_val, z)
        # return matrix also should be in RawIndex
        # since it is the order in which constraints are added

        Jac = sparse(I, J, jac_val)  # SparseMatrixCSC
        return Jac
    end

    function jacobian_callback(        
        n::Cint,
        nnz::Cint,
        z::Vector{Cdouble},
        col_start::Vector{Cint},
        col_len::Vector{Cint},
        row::Vector{Cint},
        data::Vector{Cdouble}
    )
        @assert n == length(z) == length(col_start) == length(col_len)
        @assert nnz == length(row) == length(data)    
        # z is in RawIndex, passed from PATHSolver

        Jac = j_eval(z)  # SparseMatrixCSC

        # Pouring Jac::SparseMatrixCSC to the format used in PATH
        for i in 1:n
            col_start[i] = Jac.colptr[i]
            col_len[i] = Jac.colptr[i + 1] - Jac.colptr[i]
        end

        rv = rowvals(Jac)
        nz = nonzeros(Jac)
        num_nonzeros = SparseArrays.nnz(Jac)
        for i in 1:num_nonzeros
            row[i] = rv[i]
            data[i] = nz[i]
        end

        return Cint(0)
    end

    mcp_data = get_MCP_data(m)
    n = length(mcp_data)

    # Two Indices
    # MCP_Index: the order stored in MCPModel = array index of Array{ComplementarityType}
    # RawIndex: the order used in JuMP / MathOptInteface

    # Declaring MCP mapping F as constraints
    # in order to query Jacobian using AutoDiff thru MathProgBase
    # i = RawIndex
    # Add constraint in the order of RawIndex
    p = sortperm_MCP_data(mcp_data)
    JuMP.@NLconstraint(m, [i=1:n], mcp_data[p[i]].F == 0)

    # lb and ub in RawIndex
    lb, ub = get_bounds_in_raw_index(mcp_data)

    # var_name, F_name in RawIndex
    var_name, F_name = get_names_in_raw_index(mcp_data)

    # initial values
    initial_values = get_initial_values_in_raw_index(m, mcp_data)

    # overestimating number of nonzeros in Jacobian
    nnz = max( SparseArrays.nnz(j_eval(lb)), SparseArrays.nnz(j_eval(ub)) )
    nnz = max( nnz, SparseArrays.nnz(j_eval(initial_values)) )
    for i in 1:2
        z_rand = max.(lb, min.(ub, rand(Float64, n)))
        nnz_rand = SparseArrays.nnz(j_eval(z_rand))
        nnz = max( nnz, nnz_rand )
    end
    nnz = min( 2 * nnz, n^2 )


    # Solve the MCP using PATHSolver
    status, z, info = PATHSolver.solve_mcp(
        function_callback, 
        jacobian_callback, 
        lb, 
        ub, 
        initial_values; 
        nnz = nnz, 
        variable_name = var_name,
        constraint_name = F_name,
        kwargs...
    )
    # z is in RawIndex

    # After solving set the values in m::JuMP.Model to the solution obtained.
    for i in 1:n
        set_result_value(mcp_data[i], z[mcp_data[i].raw_idx])
    end

    # Cleanup. Remove all dummy @NLconstraints added,
    # so that the model can be re-used for multiple runs
    m.nlp_data.nlconstr = []

    # This function has changed the content of m already.
    return status
end


function _solve_nlsolve!(m::JuMP.Model; method=:trust_region)
    d = JuMP.NLPEvaluator(m)

    function function_callback!(fvec, z)
        # z is in RawIndex, passed from PATHSolver
        MOI.initialize(d, [:Jac])
        F_val = zeros(n)
        MOI.eval_constraint(d, F_val, z)

        copyto!(fvec, F_val)
    end

    function jacobian_callback!(fjac, z)
        # z is in RawIndex, passed from PATHSolver
        MOI.initialize(d, [:Jac])
        J_struct = MOI.jacobian_structure(d)

        I = first.(J_struct)
        J = last.(J_struct)

        jac_val = zeros(length(J))
        MOI.eval_constraint_jacobian(d, jac_val, z)

        sparse_fjac = sparse(I, J, jac_val)
        copyto!(fjac, full(sparse_fjac))
    end

    mcp_data = get_MCP_data(m)
    n = length(mcp_data)

    # Two Indices
    # MCP_Index: the order stored in MCPModel = array index of Array{ComplementarityType}
    # RawIndex: the order used in JuMP / MathProgBase

    # Declaring MCP mapping F as constraints
    # in order to query Jacobian using AutoDiff thru MathProgBase
    # i = RawIndex
    # Add constraint in the order of RawIndex
    p = sortperm_MCP_data(mcp_data)
    JuMP.@NLconstraint(m, [i=1:n], mcp_data[p[i]].F == 0)

    # lb and ub in RawIndex
    lb, ub = get_bounds_in_raw_index(mcp_data)

    # initial values
    initial_values = get_initial_values_in_raw_index(m, mcp_data)

    # Solve the MCP using NLsolve
    # ALL inputs to NLsolve must be in RawIndex

    r = NLsolve.mcpsolve(function_callback!, jacobian_callback!, lb, ub, initial_values, method = method,
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

    # julia> fieldnames(typeof(r))
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


    # r.zero is in RawIndex
    for i in 1:n
        set_result_value(mcp_data[i], r.zero[mcp_data[i].raw_idx])
    end


    # Cleanup. Remove all dummy @NLconstraints added,
    # so that the model can be re-used for multiple runs
    m.nlp_data.nlconstr =  []

    # This function has changed the content of m already.
    return r
end


function set_result_value(mcp_data::ComplementarityType, value::Float64)
    mcp_data.result_value = value
end

function set_result_value(v::JuMP.VariableRef, value::Float64)
    mcp_data = get_MCP_data(v.model)
    for i in 1:length(mcp_data)
        if mcp_data[i].var == v
            mcp_data.result_value = value
            break
        end
    end
end

function result_value(v::JuMP.VariableRef)
    mcp_data = get_MCP_data(v.model)
    result_value = NaN
    for i in 1:length(mcp_data)
        if mcp_data[i].var == v
            result_value = mcp_data[i].result_value
            break
        end
    end
    return result_value
end



# https://stackoverflow.com/questions/50084877/how-to-alias-a-macro-in-julia#comment87277306_50085297
@eval const $(Symbol("@mapping")) = $(Symbol("@NLexpression"))


function add_complementarity(m::JuMP.Model, var::JuMP.VariableRef, F::JuMP.NonlinearExpression, F_name::String)

  lb = JuMP.has_lower_bound(var) ? JuMP.lower_bound(var) : -Inf
  ub = JuMP.has_upper_bound(var) ? JuMP.upper_bound(var) : Inf
  var_name = JuMP.name(var)
  new_dimension = ComplementarityType(lb, var, ub, F, raw_index(var), var_name, F_name, NaN)
  mcp_data = get_MCP_data(m)
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
    # if isa($F, JuMP.JuMPArray) || isa($F, Array)
    #   @assert length($F) == length($var)
    # end

    # when var is a single JuMP variable
    if isa($var, JuMP.VariableRef)
      ex_var = Meta.parse(JuMP.name($var))
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
    elseif isa($var, Array{JuMP.VariableRef})
      for idx in CartesianIndices(size($var))
        idx_name = idx.I
        F_name = string($F_base_name, "[", idx_name, "]")
        var_idx = $var[idx]
        F_idx = $F[idx]

        add_complementarity($m, var_idx, F_idx, F_name)
      end

    else # isa($var, JuMP.JuMPArray) && length(($var).indexsets) == 1 or > 1
      # when var is a multi-dimensional JuMP variable array
      ex = :(Base.product())
      for i in 1:length($var.axes)
        push!(ex.args, $var.axes[i])
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
