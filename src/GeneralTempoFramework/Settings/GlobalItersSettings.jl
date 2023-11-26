#--------------------------------------------------------------------------------------------------------------        
struct GlobalIterationsKeys
    iterative_mode::Bool
    save_global_iterations::Bool
    GlobalIterationsKeys(;iterative_mode = true, save_global_iterations = false) = new(iterative_mode, save_global_iterations)
end

# Implementation of the show method for BasicTempoKeys
function Base.show(io::IO, keys::GlobalIterationsKeys)
    println(io, "Global iterations Tempo keys:")
    print(io, "        Iterative mode: ", keys.iterative_mode)
    print(io, "        Save global iterations: ", keys.save_global_iterations)
    return nothing
end


# Настройки для глобальных итераций
struct GlobalIterationsSettings
    keys::GlobalIterationsKeys
    iters::Int64
    nits::Vector{Int64}
    gain::Vector{Float64}
    tparams_local::Vector{Vector{GeneralTempoParameter}}
end

function Base.show(io::IO, gisets::GlobalIterationsSettings)
    println(io, "Global iteration settings:")
    println(io, "   ", gisets.keys)
    println(io, "   Number of iterations: ", gisets.iters)
    for iter in 1:gisets.iters
        println(io, "   Step # $iter:")
        println(io, "       GAIN value: ", iter <= length(gisets.gain) ? gisets.gain[iter] : gisets.gain[end])
        println(io, "       NITS value: ", iter <= length(gisets.nits) ? gisets.nits[iter] : gisets.nits[end])
        println(io, "       Local Tempo parameters: ", gisets.tparams_local[iter])
    end
	return nothing
end

GlobalIterationsSettings(;
    keys = GlobalIterationsKeys(),
    iters,
    nits = 3 * ones(Int64, iters),
    gain = ones(Float64, iters),
    tparams_local = [Vector{GeneralTempoParameter}() for _ in 1:iters]
    ) = GlobalIterationsSettings(
        keys,
        iters,
        nits,
        gain,
        tparams_local
        )
