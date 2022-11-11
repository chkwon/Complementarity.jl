# isdefined(Base, :__precompile__) && __precompile__()

# https://github.com/StructJuMP/StructJuMP.jl/blob/master/src/StructJuMP.jl
# using JuMP # To reexport, should be using (not import)
# # Macro to exportall
# macro exportall(pkg)
#     Expr(:export, names(JuMP)...)
# end
# @exportall JuMP

module Complementarity

# package code goes here
# importall JuMP
using Base.Meta
using LinearAlgebra, SparseArrays

import PATHSolver, NLsolve, MathOptInterface#, OrderedCollections
const MOI = MathOptInterface

using JuMP
macro exportall(pkg)
    Expr(:export, names(JuMP)...)
end
@exportall JuMP

export  MCPModel, MCPData, ComplementarityType,
        complements, solveMCP, solveLCP,
        result_value, set_start_value,
        @complementarity, @complements, @mapping, @variable,
        @NLexpression, @expression
        PATHSolver

mpec_tol = 1e-8


include("mcp.jl")
include("mpec.jl")

end # module
