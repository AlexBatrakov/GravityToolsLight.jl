const G_CAV =  6.67408e-8
const M_sun = 1.3271244e26 / G_CAV
const c = 2.99792458e10
const d = 3600*24
const rad = 180/pi

abstract type GravityTheory end

struct GR <: GravityTheory
end

abstract type ScalarTensorGravity <: GravityTheory end

struct DEF <: ScalarTensorGravity
	alpha0::Float64
	beta0::Float64
end

DEF() = DEF(0.0, 0.0)

struct DEFGrid
    eosname::Symbol
    dims::NTuple{N,Int64} where N
    alpha0::Vector{Float64}
    beta0::Vector{Float64}
    alphaA::Array{Float64,3}
    betaA::Array{Float64,3}
    kA::Array{Float64,3}
    mA::Array{Float64,3}
end

function read_DEFGrid(eosname::Symbol, path_to_grids=".")
    path = path_to_grids * "/" * string(eosname) * "/"
    N_pc, N_alpha0, N_beta0 = readdlm(path*"dims.dat")[2,:]
    alpha0 = readdlm(path*"alpha0.dat")[:,1]
    beta0 = readdlm(path*"beta0.dat")[:,1]
    alphaA = Array{Float64}(undef, N_pc, N_alpha0, N_beta0)
    betaA = Array{Float64}(undef, N_pc, N_alpha0, N_beta0)
    kA = Array{Float64}(undef, N_pc, N_alpha0, N_beta0)
    mA = Array{Float64}(undef, N_pc, N_alpha0, N_beta0)
    alphaA_temp = readdlm(path*"alphaA.dat")
    betaA_temp = readdlm(path*"betaA.dat")
    kA_temp = readdlm(path*"kA.dat")
    mA_temp = readdlm(path*"massA.dat")
    for i_alpha0 in 1:N_alpha0
        for i_beta0 in 1:N_beta0
            alphaA[:,i_alpha0,i_beta0] .= alphaA_temp[(i_alpha0-1)*N_beta0+i_beta0,3:end]
            betaA[:,i_alpha0,i_beta0] .= betaA_temp[(i_alpha0-1)*N_beta0+i_beta0,3:end]
            kA[:,i_alpha0,i_beta0] .= kA_temp[(i_alpha0-1)*N_beta0+i_beta0,3:end]
            mA[:,i_alpha0,i_beta0] .= mA_temp[(i_alpha0-1)*N_beta0+i_beta0,3:end]
        end
    end
    return DEFGrid(eosname, (N_pc, N_alpha0, N_beta0), alpha0, beta0, alphaA, betaA, kA, mA)
end

struct DEFMassGrid
    eosname::Symbol
    N_pc::Int64
    alpha0::Float64
    beta0::Float64
    alphaA::Vector{Float64}
    betaA::Vector{Float64}
    kA::Vector{Float64}
    mA::Vector{Float64}
end

function interpolate_DEFMassGrid(grid::DEFGrid, alpha0, beta0)

    if !(grid.alpha0[1] <= alpha0 <= grid.alpha0[end]) || !(grid.beta0[1] <= beta0 <= grid.beta0[end])
        error("interpolation for alpha0=$alpha0, beta0=$beta0 is out of possible range")
    end

#    i_alpha0 = searchsortedfirst(grid.alpha0, alpha0)
#    i_beta0 = searchsortedfirst(grid.beta0, beta0)
#    i_alpha0 = i_alpha0 == grid.dims[2] ? i_alpha0 - 1 : i_alpha0
#    i_beta0 = i_beta0 == grid.dims[3] ? i_beta0 - 1 : i_beta0

    i_alpha0, i_beta0 = 1, 1
    for i in 1:grid.dims[2]-1
        if grid.alpha0[i] <= alpha0 
            i_alpha0 = i
        end
    end
    for i in 1:grid.dims[3]-1
        if grid.beta0[i] <= beta0 
            i_beta0 = i
        end
    end

    println("i_alpha0 = $i_alpha0, i_beta0 = $i_beta0")

    x_alpha0 = (alpha0 - grid.alpha0[i_alpha0]) / (grid.alpha0[i_alpha0+1] - grid.alpha0[i_alpha0])
    x_beta0 = (beta0 - grid.beta0[i_beta0]) / (grid.beta0[i_beta0+1] - grid.beta0[i_beta0])

    alphaA = grid.alphaA[:, i_alpha0, i_beta0] * (1-x_alpha0) * (1-x_beta0) + 
             grid.alphaA[:, i_alpha0+1, i_beta0] * x_alpha0 * (1-x_beta0) +
             grid.alphaA[:, i_alpha0, i_beta0+1] * (1-x_alpha0) * x_beta0 +
             grid.alphaA[:, i_alpha0+1, i_beta0+1] * x_alpha0 * x_beta0

    betaA = grid.betaA[:, i_alpha0, i_beta0] * (1-x_alpha0) * (1-x_beta0) + 
            grid.betaA[:, i_alpha0+1, i_beta0] * x_alpha0 * (1-x_beta0) +
            grid.betaA[:, i_alpha0, i_beta0+1] * (1-x_alpha0) * x_beta0 +
            grid.betaA[:, i_alpha0+1, i_beta0+1] * x_alpha0 * x_beta0

    kA = grid.kA[:, i_alpha0, i_beta0] * (1-x_alpha0) * (1-x_beta0) + 
            grid.kA[:, i_alpha0+1, i_beta0] * x_alpha0 * (1-x_beta0) +
            grid.kA[:, i_alpha0, i_beta0+1] * (1-x_alpha0) * x_beta0 +
            grid.kA[:, i_alpha0+1, i_beta0+1] * x_alpha0 * x_beta0

    mA = grid.mA[:, i_alpha0, i_beta0] * (1-x_alpha0) * (1-x_beta0) + 
            grid.mA[:, i_alpha0+1, i_beta0] * x_alpha0 * (1-x_beta0) +
            grid.mA[:, i_alpha0, i_beta0+1] * (1-x_alpha0) * x_beta0 +
            grid.mA[:, i_alpha0+1, i_beta0+1] * x_alpha0 * x_beta0

    return DEFMassGrid(grid.eosname, grid.dims[1], alpha0, beta0, alphaA, betaA, kA, mA)
end

function interpolate_NS(mgrid::DEFMassGrid, mA)
    if !(mgrid.mA[1] <= mA <= mgrid.mA[end])
        error("interpolation for mA=$mA is out of possible range")
    end

    i_mA = 1
    for i in 1:mgrid.N_pc
        if mgrid.mA[i] <= mA
            i_mA = i
        end
    end

    x_mA = (mA - mgrid.mA[i_mA]) / (mgrid.mA[i_mA+1] - mgrid.mA[i_mA])
    alphaA = mgrid.alphaA[i_mA] * (1-x_mA) + mgrid.alphaA[i_mA+1] * x_mA
    betaA = mgrid.betaA[i_mA] * (1-x_mA) + mgrid.betaA[i_mA+1] * x_mA
    kA = mgrid.kA[i_mA] * (1-x_mA) + mgrid.kA[i_mA+1] * x_mA
    return (alphaA = alphaA, betaA = betaA, kA = kA)
end

KType = NamedTuple{(:Pb, :T0, :e0, :omega0, :x0), NTuple{5, Float64}}

PKType = NamedTuple{(:k, :gamma, :Pbdot, :r, :s), NTuple{5, Float64}}

mutable struct Object
    type::Symbol
    mass::Float64
    alphaA::Float64
    betaA::Float64
    kA::Float64
    function Object(type::Symbol)
        return new(type)
    end
    function Object(type::Symbol, mass::Float64)
        return new(type, mass)
    end
end

mutable struct BinarySystem
    PSR::Object
    Comp::Object
    K_params::KType
    PK_params::PKType
    function BinarySystem(PSR::Object, Comp::Object)
        return new(PSR, Comp)
#        return new(PSR, Comp, KType(zeros(5)), PKType(zeros(5)))
    end
    function BinarySystem(PSR_type::Symbol, Comp_type::Symbol)
#        return new(Object(PSR_type), Object(Comp_type), KType(zeros(5)), PKType(zeros(5)))
        return new(Object(PSR_type), Object(Comp_type))
    end
end

function Base.show(io::IO, bnsys::BinarySystem)
	println(io, "Pulsar: ", bnsys.PSR)
	println(io, "Companion: ", bnsys.Comp)
	println(io, "Keplerian parameters: ", bnsys.K_params)
	println(io, "Post-Keplerian parameters: ", bnsys.PK_params)
	return nothing
end

struct Settings
    path_to_grids::String
end

abstract type PhysicalFramework end

mutable struct GRPhysicalFramework <: PhysicalFramework
    theory::GR
    bnsys::BinarySystem
end

mutable struct DEFPhysicalFramework <: PhysicalFramework
    theory::DEF
    eosname::Symbol
    grid::DEFGrid
    bnsys::BinarySystem
    sets::Settings
    mgrid::DEFMassGrid
    function DEFPhysicalFramework(theory::DEF, eosname::Symbol, bnsys::BinarySystem, sets::Settings)
        grid = read_DEFGrid(eosname, sets.path_to_grids)
        return(new(theory, eosname, grid, bnsys, sets))
    end
end

function Base.show(io::IO, pf::DEFPhysicalFramework)
	println(io, "DEFPhysicalFramework:")
	println(io, "Theory: ", pf.theory)
	println(io, "EOS name: ", pf.eosname)
	println(io, "Binary system:")
	println(io, pf.bnsys)
	return nothing
end

function reinitialize!(pf::DEFPhysicalFramework)
    pf.grid = read_DEFGrid(pf.eosname, pf.sets.path_to_grids)
    return pf
end

#function calculate!(bnsys::BinarySystem)
#end