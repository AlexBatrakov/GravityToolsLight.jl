struct AnnealingParameter
    name::String
    initial_value::Float64
    min_value::Float64
    max_value::Float64
    initial_step_size::Float64
    is_angle::Bool
end

const AP = AnnealingParameter

struct AnnealingSettings
    parameters::Vector{AnnealingParameter}
    initial_temperature::Float64
    cooling_rate::Float64
    minimum_temperature::Float64
    iterations_per_temp::Int
    max_iterations::Int
    parallel::Bool
    random_seed::Int
end

AnnealingSettings(;parameters, initial_temperature, cooling_rate, minimum_temperature, iterations_per_temp, max_iterations, parallel) = 
    AnnealingSettings(parameters, initial_temperature, cooling_rate, minimum_temperature, iterations_per_temp, max_iterations, parallel, rand(1:10000))


