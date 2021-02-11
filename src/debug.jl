

# package code goes here
# importall JuMP
using Base.Meta
using LinearAlgebra, SparseArrays

import PATHSolver, NLsolve, MathOptInterface
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
        @NLexpression, @expression,
        PATHSolver

mpec_tol = 1e-8


include("mcp.jl")
include("mpec.jl")
