#--------------------------------------------------------------------------------------------------------------        
# Настройки для глобальных итераций
struct GlobalIterationsSettings
    iters::Int64
    nits::Vector{Int64}
    gain::Vector{Float64}
    tparams_local::Vector{Vector{GeneralTempoParameter}}
end

function Base.show(io::IO, gisets::GlobalIterationsSettings)
    println(io, "Global iteration settings:")
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
    iters,
    nits = 3 * ones(Int64, iters),
    gain = ones(Float64, iters),
    tparams_local = [Vector{GeneralTempoParameter}() for _ in 1:iters]
    ) = GlobalIterationsSettings(
        iters,
        nits,
        gain,
        tparams_local
        )
