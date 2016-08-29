# isdefined(Base, :__precompile__) && __precompile__()

module Complementarity

# package code goes here

using JuMP
using Base.Meta

import PATHSolver, MathProgBase, NLsolve
importall JuMP, Ipopt

export  MCPModel, MCPData, Complementarity,
        complements, solveMCP,
        @complements

include("mcp.jl")
include("mpec.jl")

end # module
