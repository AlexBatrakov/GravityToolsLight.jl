mutable struct PriorSettings
    prior_name::String
    prior_file::String
    prior_values::Vector{Float64}
    prior_quntile_targets::Vector{Float64}
end

function PriorSettings(;prior_name::String, prior_file::String, prior_quntile_targets::Vector{Float64})
    prior_values = readdlm(prior_file)[1:end]
    return PriorSettings(prior_name, prior_file, prior_values, prior_quntile_targets)
end

