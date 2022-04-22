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

struct DEFSmallGrid
    eosname::Symbol
    N_pc::Int64
    alpha0::Float64
    beta0::Float64
    alphaA::Vector{Float64}
    betaA::Vector{Float64}
    kA::Vector{Float64}
    mA::Vector{Float64}
end

function interpolate_DEFSmallGrid(grid::DEFGrid, alpha0, beta0)

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

    return DEFSmallGrid(grid.eosname, grid.dims[1], alpha0, beta0, alphaA, betaA, kA, mA)
end

function interpolate_in_mass(sgrid::DEFSmallGrid, mA)
    if !(sgrid.mA[1] <= mA <= sgrid.mA[end])
        error("interpolation for mA=$mA is out of possible range")
    end

    i_mA = 1
    for i in 1:sgrid.N_pc
        if sgrid.mA[i] <= mA
            i_mA = i
        end
    end

    x_mA = (mA - sgrid.mA[i_mA]) / (sgrid.mA[i_mA+1] - sgrid.mA[i_mA])
    alphaA = sgrid.alphaA[i_mA] * (1-x_mA) + sgrid.alphaA[i_mA+1] * x_mA
    betaA = sgrid.betaA[i_mA] * (1-x_mA) + sgrid.betaA[i_mA+1] * x_mA
    kA = sgrid.kA[i_mA] * (1-x_mA) + sgrid.kA[i_mA+1] * x_mA
    return alphaA, betaA, kA
end