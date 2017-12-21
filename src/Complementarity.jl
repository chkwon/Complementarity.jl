# isdefined(Base, :__precompile__) && __precompile__()
using JuMP


module Complementarity

# package code goes here
# importall JuMP
using JuMP
using Base.Meta

import PATHSolver, MathProgBase, NLsolve

export  MCPModel, MCPData, ComplementarityType,
        complements, solveMCP, solveLCP,
        @complementarity, @complements, @mapping

mpec_tol = 1e-8


include("mcp.jl")
include("mpec.jl")

end # module
