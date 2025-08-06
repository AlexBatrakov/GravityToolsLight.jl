struct AnnealingParameter
    name::String
    initial_value::Float64
    min_value::Float64
    max_value::Float64
    initial_step_size::Float64
    is_angle::Bool
end

const AP = AnnealingParameter

function update_parameter_value(param::AnnealingParameter, new_value::Float64)
    if param.is_angle
        # Приводим значение угла к диапазону [0, param.max_value)
        new_value = mod(new_value, param.max_value)
    else
        # Это обычный параметр, не угол
        new_value = clamp(new_value, param.min_value, param.max_value)
    end
    return new_value
end

struct AnnealingSettings
    parameters::Vector{AnnealingParameter}
    energy_scale::Float64
    quenching_factor::Float64
    initial_temperature::Float64
    cooling_rate::Float64
    minimum_temperature::Float64
    step_law::Float64
    iterations_per_temp::Int
    max_iterations::Int
    parallel::Bool
    gradient::Bool
    random_seed::Int
end

AnnealingSettings(;parameters, energy_scale, quenching_factor, initial_temperature, cooling_rate, minimum_temperature, step_law, iterations_per_temp, max_iterations, parallel, gradient) = 
    AnnealingSettings(parameters, energy_scale, quenching_factor, initial_temperature, cooling_rate, minimum_temperature, step_law, iterations_per_temp, max_iterations, parallel, gradient, rand(1:10000))


