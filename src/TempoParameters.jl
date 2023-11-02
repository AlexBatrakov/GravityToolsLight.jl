#--------------------------------------------------------------------------------------------------------------
# Setting for tempo parameters

# Set precision for BigFloat to 80 digits
setprecision(BigFloat, 80)

# Definition of a general tempo parameter structure
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

# Constructor that parses a string line to extract GeneralTempoParameter
function extract_GeneralTempoParameter(line::String)
    # логика тут такая. Кажждая строчка выглядит так
    # имя_параметра значение флаг погрешность

    line_split = split(line)
    n = length(line_split)
    line_parsed = parse_tparam_field.(line_split)
    line_parsed_types = typeof.(line_parsed)


    # имя параметра может быть из одного слова или из трех
    if n >= 3 && line_parsed_types[1:3] == [String, String, String]
        n_name = 3
        name = join(line_split[1:3], " ")
    else
        n_name = 1
        name =  String(line_split[1])
    end

    # значение всегда стоит после имени долину которого мы уже определили
    value = n_name < n ? line_parsed[n_name + 1] : nothing

    # после значения обычно стоит флаг, только если следующее число не является BigFloat
    flag = n_name + 1 < n && !isa(line_parsed[n_name + 2], BigFloat) ? line_parsed[n_name + 2] : nothing

    #после флага или сразу после значения может стоять погрешность, если она отсуствует, то она не определена
    uncertainty = n_name + 1 < n && isa(line_parsed[n_name + 2], BigFloat) ? line_parsed[n_name + 2] : n_name + 2 < n ? line_parsed[n_name + 3] : nothing

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

# function format_for_value(value::Number)
#     if abs(value) < 1e-3 || abs(value) > 1e3
#         return "%.16e"
#     else
#         return "%.16f"
#     end
# end

# format_for_value(value::String) = "%s"

# function print_tparam(tp::GeneralTempoParameter)
#     val_fmt = format_for_value(tp.value)
#     unc_fmt = format_for_value(tp.uncertainty)
    
#     val_str = @sprintf(val_fmt, tp.value)
#     unc_str = tp.uncertainty !== nothing ? @sprintf(unc_fmt, tp.uncertainty) : ""

#     println("$(rpad(tp.name, 15)) $(lpad(val_str, 25)) $(tp.flag) $(lpad(unc_str, 25))")
# end


#--------------------------------------------------------------------------------------------------------------
# tempo par files

mutable struct TempoParFile
    name::String
    tparams::Dict{Symbol,GeneralTempoParameter}
    order::Vector{Symbol}
end

function TempoParFile(par_file_path::String; new_name_suffix="")
    par_file = TempoParFile(par_file_path, Dict{Symbol,GeneralTempoParameter}(), Vector{Symbol}())
    read_par_file!(par_file)
    par_file.name = par_file.name[1:end-4] * new_name_suffix * ".par"
    return par_file
end

function Base.show(io::IO, par_file::TempoParFile)
    println(io, "Tempo parameter file $(par_file.name): ")
    for (i, name_symbol) in enumerate(par_file.order)
        #print(IOContext(io, :indent => indent+4), par_file.tparams[name_symbol])
        print("    ", get_par_file_representation(par_file.tparams[name_symbol]))
        if i < length(par_file.order)
            print("\n")
        end
    end
	return nothing
end

function read_par_file!(par_file::TempoParFile)
    par_file.order = Vector{Symbol}()
    open(par_file.name, "r") do file_in
        for line in eachline(file_in)
            if startswith(line, "C ") || startswith(line, "c ")
                continue
            end
            tparam = extract_GeneralTempoParameter(line)
            par_file.tparams[tparam.name_symbol] = tparam
            push!(par_file.order, tparam.name_symbol)
        end
    end
    return par_file
end

function write_par_file(par_file::TempoParFile, name_out=par_file.name)
    open(name_out, "w") do file_out
        for name_symbol in par_file.order
            tparam = par_file.tparams[name_symbol]
            println(file_out, get_par_file_representation(tparam))
        end
    end
    return par_file
end

function modify_par_file!(par_file::TempoParFile, tparam_name::Symbol, value=tparam.value, flag=tparam.flag, uncertainty=tparam.uncertainty)
    return par_file
end

function update_par_file()

end

function extend_par_file!(par_file::TempoParFile, tparam::GeneralTempoParameter)
    if tparam.value == nothing
        par_file.tparams[tparam.name_symbol].flag = tparam.flag
    else
        if !haskey(par_file.tparams, tparam.name_symbol)
            push!(par_file.order, tparam.name_symbol)
        end
        par_file.tparams[tparam.name_symbol] = tparam
    end
    return par_file
end

function generate_par_file_path(base_name::String, suffix::String, work_dir::String)
    # Убедимся, что имя базы и суффикс не начинаются с точки и не содержат расширения файла
    base_name = strip_extension(base_name)
    suffix = lstrip(suffix, '.')

    # Сгенерируем полный путь к файлу
    if suffix != ""
        return joinpath(work_dir, base_name * "_" * suffix * ".par")
    else
        return joinpath(work_dir, base_name * ".par")
    end

    return par_file_path
end

# Убираем расширение файла, если оно есть
function strip_extension(file_name::String)
    return splitext(file_name)[1]
end
