# isdefined(Base, :__precompile__) && __precompile__()
using JuMP


module Complementarity

# package code goes here
# importall JuMP
using JuMP
using Base.Meta
using LinearAlgebra, SparseArrays

import PATHSolver, NLsolve, MathProgBase

export  MCPModel, MCPData, ComplementarityType,
        complements, solveMCP, solveLCP,
        getvalue, setvalue,
        @complementarity, @complements, @mapping, @variable,
        @NLexpression, @expression,
        PATHSolver

mpec_tol = 1e-8


include("mcp.jl")
include("mpec.jl")

end # module
