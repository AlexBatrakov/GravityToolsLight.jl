#--------------------------------------------------------------------------------------------------------------
# Setting for tempo parameters

# Set precision for BigFloat to 80 digits
setprecision(BigFloat, 80)

"""
    mutable struct GeneralTempoParameter

Represents a general parameter used in Tempo software.

# Fields
- `name::String`: The name of the parameter.
- `name_symbol::Symbol`: The symbol representation of the parameter name.
- `value::Union{Int64, BigFloat, String, Nothing}`: The value of the parameter, supporting various types.
- `flag::Union{Int64, Nothing}`: An optional flag associated with the parameter.
- `uncertainty::Union{BigFloat, Nothing}`: The uncertainty associated with the parameter value, if applicable.

# Constructors
- `GeneralTempoParameter(name::String, value=nothing; flag=nothing, uncertainty=nothing)`: Creates a new parameter, converting `Float64` to `BigFloat` and setting default values for flag and uncertainty.
- `GeneralTempoParameter(var::ValueVariable)`: Converts from `ValueVariable` type to `GeneralTempoParameter`.

# Methods
- `setproperty!(x::GeneralTempoParameter, f::Symbol, v::Float64)`: Assigns a `Float64` value to the parameter, converting it to `BigFloat`.
- `extract_GeneralTempoParameter(line::String)`: Parses a string line to extract a `GeneralTempoParameter`.
- `get_par_file_representation(tparam::GeneralTempoParameter)`: Returns a formatted representation of a `GeneralTempoParameter` for a `.par` file.

"""
mutable struct GeneralTempoParameter
    name::String
    name_symbol::Symbol
    value::Union{Int64, BigFloat, String, Nothing}
    flag::Union{Int64, Nothing}
    uncertainty::Union{BigFloat, Nothing}
end

TP = GeneralTempoParameter

# Constructor for converting Float64 to BigFloat and providing default arguments for flag and uncertainty
function GeneralTempoParameter(name::String, value::Union{Int64, Float64, BigFloat, String, Nothing}=nothing; flag::Union{Int64, Nothing}=nothing, uncertainty::Union{Float64, BigFloat, Nothing}=nothing)
    name_symbol = Symbol(name)  # Convert name to symbol
    big_value = isa(value, Float64) ? BigFloat(value) : value
    big_uncertainty = isnothing(uncertainty) ? nothing : BigFloat(uncertainty)
    return GeneralTempoParameter(name, name_symbol, big_value, flag, big_uncertainty)
end

function Base.setproperty!(x::GeneralTempoParameter, f::Symbol, v::Float64)
    if f === :value
        v = BigFloat(v)
    end
    setfield!(x, f, v)
end


# Constructor to convert from ValueVariable type
GeneralTempoParameter(var::ValueVariable) = GeneralTempoParameter(var.name, var.value)

# Display function for GeneralTempoParameter
function Base.show(io::IO, tparam::GeneralTempoParameter)
    indent = get(io, :indent, 0)
    print(io, " "^indent, tparam.name)

    # Formatted output for BigFloat values
    if tparam.value isa BigFloat
        print(io, " ", @sprintf("%.10f", tparam.value))
    else
        print(io, " ", tparam.value)
    end

    if tparam.flag !== nothing
        print(io, " ", tparam.flag)
    end

    # Formatted output for BigFloat uncertainties
    if tparam.uncertainty !== nothing
        if tparam.uncertainty isa BigFloat
            print(io, " ±", @sprintf("%.10f", tparam.uncertainty))
        else
            print(io, " ±", tparam.uncertainty)
        end
    end

    return nothing
end

# Function to parse a tempo parameter field
function parse_tparam_field(value_str)
    value_int64 = tryparse(Int64, value_str)
    if value_int64 !== nothing
        return value_int64
    end

    value_bigfloat = tryparse(BigFloat, value_str)
    if value_bigfloat !== nothing
        return value_bigfloat
    end
    
    return String(value_str)
end

# Parses a string line to create a GeneralTempoParameter object.
# Each line in the input string should be formatted as follows:
# parameter_name value flag uncertainty
function extract_GeneralTempoParameter(line::String)
    # Split the line into words to process each part
    line_split = split(line)
    n = length(line_split)
    line_parsed = parse_tparam_field.(line_split)
    line_parsed_types = typeof.(line_parsed)

    # The parameter name can be one word or a combination of three words
    if n >= 3 && line_parsed_types[1:3] == [String, String, String]
        n_name = 3
        name = join(line_split[1:3], " ")
    else
        n_name = 1
        name = String(line_split[1])
    end

    # The value always follows the name whose length we've already determined
    value = n_name < n ? line_parsed[n_name + 1] : nothing

    # After the value, there's usually a flag, only if the next number is not a BigFloat
    flag = n_name + 1 < n && !isa(line_parsed[n_name + 2], BigFloat) ? line_parsed[n_name + 2] : nothing

    # After the flag or right after the value, there may be an uncertainty, if it is absent, then it is undefined
    uncertainty = n_name + 1 < n && isa(line_parsed[n_name + 2], BigFloat) ? line_parsed[n_name + 2] : n_name + 2 < n ? line_parsed[n_name + 3] : nothing

    # Return a new GeneralTempoParameter object with the extracted values
    return GeneralTempoParameter(name, value, flag=flag, uncertainty=uncertainty)
end


# Utility function to align a string to a specified width
function align_str(s::String, n::Int)
    return s * " "^(n - length(s))
end

# Function to get a formatted representation of a GeneralTempoParameter for a .par file
function get_par_file_representation(tparam::GeneralTempoParameter)
    n_name = (tparam.value isa BigFloat && tparam.value < 0) ? 19 : 20
    n_value = (tparam.value isa BigFloat && tparam.value < 0) ? 30 : 29
    n_flag = 6
    n_uncertainty = 29

    line = align_str(tparam.name, n_name)


    # If the value is a BigFloat, use @sprintf for formatting
    if tparam.value isa BigFloat
        if 0.0 < abs(tparam.value) < 1e-3
            line *= align_str(@sprintf("%.21e", tparam.value), n_value)
        else 
            line *= align_str(@sprintf("%.21g", tparam.value), n_value)
        end
    else
        line *= align_str(string(tparam.value), n_value)
    end

    if tparam.flag !== nothing
        #line *= align_str(string(tparam.flag), n_flag)
        line *= string(tparam.flag) * "  "
    else
        line *= "   "
    end

    # If uncertainty is a BigFloat, use @sprintf for formatting
    if tparam.uncertainty !== nothing
        if tparam.value isa BigFloat && (0.0 < abs(tparam.value) < 1e-3)
            line *= align_str(@sprintf("%.21e", tparam.uncertainty), n_uncertainty)
        else
            line *= align_str(@sprintf("%.21f", tparam.uncertainty), n_uncertainty)
        end
    end

    return line
end

# Функция для обновления или добавления параметра в список tparams
function update_or_add_tparam(tparams::Vector{GeneralTempoParameter}, tparam::GeneralTempoParameter)
    # Поиск параметра в списке
    found = false
    for (i, existing_tparam) in enumerate(tparams)
        if existing_tparam.name == tparam.name
            tparams[i] = tparam # Обновление значения параметра
            found = true
            break
        end
    end

    # Добавление параметра, если он не был найден
    if !found
        push!(tparams, tparam)
    end

    return tparams
end

