const G_CAV =  6.67408e-8
const M_sun = 1.3271244e26 / G_CAV
const c = 2.99792458e10
const d = 3600*24
const rad = 180/pi

abstract type GravityTheory end

struct GR <: GravityTheory
end

abstract type ScalarTensorGravity <: GravityTheory end

mutable struct DEF <: ScalarTensorGravity
	alpha0::Float64
	beta0::Float64
    function DEF(alpha0, beta0)
        return new(alpha0, beta0)
    end
    function DEF()
        return new()
    end
end

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
    mAmax::Float64
    i_mAmax::Int64
    function DEFMassGrid(eosname, N_pc, alpha0, beta0, alphaA, betaA, kA, mA)
        mAmax, i_mAmax = findmax(mA)
        return new(eosname, N_pc, alpha0, beta0, alphaA, betaA, kA, mA, mAmax, i_mAmax)
    end
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

#    println("i_alpha0 = $i_alpha0, i_beta0 = $i_beta0")

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



#    if !(mgrid.mA[1] <= mA <= mgrid.mA[end])
#        error("interpolation for mA=$mA is out of possible range [$(mgrid.mA[1]), $(mgrid.mA[end])]")
#    end

    i_mA = 1
    for i in 1:mgrid.i_mAmax-1
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

K_list = (:Pb, :T0, :e0, :omega0, :x0)
PK_list = (:k, :gamma, :Pbdot, :r, :s, :h3, :varsigma, :dtheta)
X_list = (:m2, :q, :deltaN)

KType = NamedTuple{K_list, NTuple{length(K_list), Float64}}
PKType = NamedTuple{PK_list, NTuple{length(PK_list), Float64}}
XType = NamedTuple{X_list, NTuple{length(X_list), Float64}}

mutable struct Object
    type::Symbol
    mass::Float64
    alphaA::Float64
    betaA::Float64
    kA::Float64
    function Object(type::Symbol)
        return new(type, 0.0, 0.0, 0.0, 0.0)
    end
    function Object(type::Symbol, mass::Float64)
        return new(type, mass, 0.0, 0.0, 0.0)
    end
end

function Base.show(io::IO, obj::Object)
	println(io, "Type of object: ", obj.type)
	println(io, "   Mass (M⊙): ", obj.mass)
	println(io, "   alphaA: ", obj.alphaA)
	println(io, "   betaA: ", obj.betaA)
    println(io, "   kA: ", obj.kA)
	return nothing
end

mutable struct BinarySystem
    psr::Object
    comp::Object
    name::String
    K_params::KType
    PK_params::PKType
    X_params::XType
    function BinarySystem(psr::Object, comp::Object)
        return new(psr, comp)
    end
    function BinarySystem(psr_type::Symbol, comp_type::Symbol)
        return new(Object(psr_type), Object(comp_type))
    end
end

function Base.show(io::IO, bnsys::BinarySystem)
    println(io, "Binary system:")
	print(io,   "Pulsar:\n   ", bnsys.psr)
	print(io,   "Companion:\n   ", bnsys.comp)
	println(io, "Keplerian parameters:\n   ", bnsys.K_params)
	println(io, "Post-Keplerian parameters:\n   ", bnsys.PK_params)
    print(io,   "Extra parameters:\n   ", bnsys.X_params)
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
    bnsys::BinarySystem
    sets::Settings
    grid::DEFGrid
    mgrid::DEFMassGrid
    function DEFPhysicalFramework(theory::DEF, eosname::Symbol, bnsys::BinarySystem)
        return new(theory, eosname, bnsys) 
    end
    function DEFPhysicalFramework(theory::DEF, eosname::Symbol, bnsys::BinarySystem, sets::Settings)
        return new(theory, eosname, bnsys, sets)
    end
end

function Base.show(io::IO, pf::DEFPhysicalFramework)
	println(io, "DEFPhysicalFramework:")
	println(io, "Theory: ", pf.theory)
	println(io, "EOS name: ", pf.eosname)
	print(io, pf.bnsys)
	return nothing
end

function read_grid!(pf::DEFPhysicalFramework)
    pf.grid = read_DEFGrid(pf.eosname, pf.sets.path_to_grids)
    return pf
end

function interpolate_mgrid!(pf::DEFPhysicalFramework)
    if pf.theory.alpha0 == 0 && pf.theory.beta0 == 0
        println("No interpolation for the GR")
        return pf
    end
    pf.mgrid = interpolate_DEFMassGrid(pf.grid, pf.theory.alpha0, pf.theory.beta0)
    return pf
end

function interpolate_psr!(pf::DEFPhysicalFramework)
    psr = pf.bnsys.psr
    if pf.theory.alpha0 == 0.0 && pf.theory.beta0 == 0.0
        psr.alphaA, psr.betaA, psr.kA = 0.0, 0.0, 0.0
    else
        psr.alphaA, psr.betaA, psr.kA = interpolate_NS(pf.mgrid, psr.mass)
    end
    return pf
end

function interpolate_comp!(pf::DEFPhysicalFramework)
    comp = pf.bnsys.comp
    if comp.type == :NS
        if pf.theory.alpha0 == 0.0 && pf.theory.beta0 == 0.0
            comp.alphaA, comp.betaA, comp.kA = 0.0, 0.0, 0.0
        else
            comp.alphaA, comp.betaA, comp.kA = interpolate_NS(pf.mgrid, comp.mass)
        end
    elseif  comp.type == :BH
        comp.alphaA, comp.betaA, comp.kA = 0.0, 0.0, 0.0
    elseif  comp.type == :WD
        comp.alphaA, comp.betaA, comp.kA = pf.theory.alpha0, pf.theory.beta0, 0.0
    else
        error("the type $(comp.type) of the compnion is not supported")
    end
    return pf
end

function calculate_PK_params!(pf::DEFPhysicalFramework)
	alpha0 = pf.theory.alpha0
	beta0 = pf.theory.beta0
	Pb = pf.bnsys.K_params.Pb * d
	e = pf.bnsys.K_params.e0
	x = pf.bnsys.K_params.x0
	m1 = pf.bnsys.psr.mass * M_sun
	m2 = pf.bnsys.comp.mass * M_sun
	m = m1 + m2

	m2_bare = m2/(1+alpha0^2)
	alphaA, betaA, kA = pf.bnsys.psr.alphaA, pf.bnsys.psr.betaA, pf.bnsys.psr.kA
    alphaB, betaB, kB = pf.bnsys.comp.alphaA, pf.bnsys.comp.betaA, pf.bnsys.comp.kA

	G = G_CAV / (1+alpha0^2)

	GAB = G*(1 + alphaA*alphaB)

    GAB = abs(GAB) #the check of the sign

	n = 2*pi/Pb

    gamma=(e*m2*(m + alphaB*kA*m + m2 + alphaA*alphaB*m2)*((GAB*m*n)/c^3)^(2/3))/(m^2*(n + alphaA*alphaB*n))  

    k = (((-6 + alphaA*(2*alphaB*(-2 + alphaA*alphaB) + alphaA*betaB))*m + (alphaB^2*betaA - alphaA^2*betaB)*m2)*((GAB*m*n)/c^3)^(2/3))/(2. *(1 + alphaA*alphaB)^2*(-1 + e^2)*m)


    s = (c*m*x)/(((GAB*m)/n^2)^(1/3)*(m2 - (m2*(18 - (20*alphaA*alphaB)/(1 + alphaA*alphaB) - (2*m2)/m + (2*m2^2)/m^2)*((GAB*m*n)/c^3)^(2/3))/6.))

#    s = c*(m*n*x)/(m2*(GAB*m*n)^(1/3))

    Pbdot_m = (-3*e^2*(4 + e^2)*m1*m2*(alphaA + alphaB*((5/3) + betaA/(1 + alphaA*alphaB)) + (alphaA*betaB)/(1 + alphaA*alphaB) + (2*(alphaA - alphaB)*m2)/(3. *m))^2*((GAB*m*n)/c^3)^(5/3)*pi)/(4. *(1 + alphaA*alphaB)*(1 - e^2)^3.5*m^2)

    Pbdot_d = (2*m1*m2*(((alphaA - alphaB)^2*(-2 + e^2 + e^4)*GAB*m*n)/(2. *c^3) - (2*(alphaA - alphaB)*(((1 + alphaA*alphaB)*(32 + 124*e^2 + 19*e^4)*(m - 2*m2)*(alphaA*m1 + alphaB*m2))/ 4. + 5*(1 + 3*e^2 + (3*e^4)/8.)*m*(alphaA*betaB*m1 - alphaB*betaA*m2))* ((GAB*m*n)/c^3)^(5/3))/(5. *(1 + alphaA*alphaB)*m^2))*pi)/ ((1 + alphaA*alphaB)*(1 - e^2)^3.5*m^2)

    Pbdot_qg = (-2*(96 + 292*e^2 + 37*e^4)*m1*m2*((GAB*m*n)/c^3)^(5/3)*pi)/(5. *(1 - e^2)^3.5*m*(m + alphaA*alphaB*m))

    Pbdot_qphi = -((96 + 292*e^2 + 37*e^4)*m1*m2*(alphaB*m1 + alphaA*m2)^2*((GAB*m*n)/c^3)^(5/3)*pi)/(15. *(1 + alphaA*alphaB)*(1 - e^2)^3.5*m^4)

    Pbdot = Pbdot_m + Pbdot_d + Pbdot_qg + Pbdot_qphi

#    println(alphaA, " ", betaA, " ", kA)
#    println(alphaB, " ", betaB, " ", kB)
#    println(m1/M_sun, " ", m2/M_sun, " ", Pbdot_m, " ", Pbdot_d, " ", Pbdot_qg, " ", Pbdot_qphi)

#    r =  G*m2_bare / c^3
	r =  G*m2 / c^3

    varsigma = (1.0 - sqrt(1 - s^2 >= 0 ? 1 - s^2 : NaN))/s
    h3 = r * varsigma^3

    dtheta = (GAB*m*n/c^3)^(2/3) / (m^2 * (1+alphaA*alphaB)) * ((3.5-0.5*alphaA*alphaB)*m1^2 + (6-alphaA*alphaB-kA*alphaB)*m1*m2 + (2-alphaA*alphaB-kA*alphaB)*m2^2)

    pf.bnsys.PK_params = (k = k, gamma = gamma, Pbdot = Pbdot, r = r, s = s, h3 = h3, varsigma = varsigma, dtheta = dtheta)

	return pf
end

function calculate_X_params!(pf::DEFPhysicalFramework)
    m2 = pf.bnsys.comp.mass
    q = pf.bnsys.psr.mass / pf.bnsys.comp.mass
    deltaN = pf.theory.alpha0*(pf.bnsys.psr.alphaA - pf.theory.alpha0)
    pf.bnsys.X_params = (m2 = m2, q = q, deltaN = deltaN)
	return pf
end



function interpolate_bnsys!(pf::DEFPhysicalFramework)
    interpolate_psr!(pf)
    interpolate_comp!(pf)
    calculate_PK_params!(pf)
    calculate_X_params!(pf)
    return pf
end

#function calculate!(bnsys::BinarySystem)
#end