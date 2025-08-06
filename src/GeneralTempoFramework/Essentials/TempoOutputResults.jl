#--------------------------------------------------------------------------------------------------------------

abstract type AbstractTempoResult end

mutable struct BasicTempoOutputResult <: AbstractTempoResult
    rms_pre_fit_residual_us::Float64
    rms_post_fit_residual_us::Float64
    rms_tn_post_fit_residual_us::Float64
    chisqr::Float64
    nfree::Int
    chisqr_red::Float64
    pre_post::Float64
    number_of_fit_parameters::Int
    number_of_points_in_fit::Int
    offset::Tuple{Float64, Float64}
    offset_e_sqrt_n::Float64
    AD_objective::Float64
    chisqr_full::Float64
end

BasicTempoOutputResult() = BasicTempoOutputResult(NaN, NaN, NaN, NaN, 0, NaN, NaN, 0, 0, (NaN, NaN), NaN, NaN, NaN)

function Base.show(io::IO, result::BasicTempoOutputResult)
    indent = get(io, :indent, 0)
    println(io, ' '^indent, "Basic Tempo Output Result:")
    println(io, ' '^(indent+4), "RMS Pre-fit Residual (us): ", result.rms_pre_fit_residual_us)
    println(io, ' '^(indent+4), "RMS Post-fit Residual (us): ", result.rms_post_fit_residual_us)
    println(io, ' '^(indent+4), "RMS TN Post-fit Residual (us): ", result.rms_tn_post_fit_residual_us)
    println(io, ' '^(indent+4), "Chisqr: ", result.chisqr)
    println(io, ' '^(indent+4), "Nfree: ", result.nfree)
    println(io, ' '^(indent+4), "Chisqr/Nfree: ", result.chisqr_red)
    println(io, ' '^(indent+4), "Pre/Post: ", result.pre_post)
    println(io, ' '^(indent+4), "Number of Fit Parameters: ", result.number_of_fit_parameters)
    println(io, ' '^(indent+4), "Number of Points in Fit: ", result.number_of_points_in_fit)
    println(io, ' '^(indent+4), "Offset: (", result.offset[1], ", ", result.offset[2], ")")
    println(io, ' '^(indent+4), "Offset E * Sqrt(n): ", result.offset_e_sqrt_n)
    println(io, ' '^(indent+4), "AD objective: ", result.AD_objective)
    println(io, ' '^(indent+4), "Chisqr_full: ", result.chisqr_full)
end

struct FitParameter
    name::String
    name_symbol::Symbol
    pre_fit::Float64
    post_fit::Float64
    uncertainty::Float64
    difference::Float64
    fit_flag::Bool

    function FitParameter(name::String, pre_fit::Float64, post_fit::Float64, uncertainty::Float64, difference::Float64, fit_flag::Bool)
        new(name, Symbol(name), pre_fit, post_fit, uncertainty, difference, fit_flag)
    end
end


function Base.show(io::IO, param::FitParameter)
    indent = get(io, :indent, 0)
    println(io, ' '^(indent), param.name, "  Pre-fit: ", param.pre_fit, "  Post-fit: ", param.post_fit, "  Unc: ", param.uncertainty, "  Diff: ", param.difference, "  fit: ", param.fit_flag)
    # и так далее для остальных полей
end

struct DetailedTempoOutputResult <: AbstractTempoResult
    basic::BasicTempoOutputResult
    fit_parameters::StructArray{FitParameter}
    fit_parameters_order::Dict{Symbol, Int64}

    function DetailedTempoOutputResult(basic::BasicTempoOutputResult, fit_parameters::Vector{FitParameter})
        # Создание словаря индексов
        order = Dict(Symbol(param.name) => i for (i, param) in enumerate(fit_parameters))

        # Создание StructArray из вектора параметров
        params_struct_array = StructArray(fit_parameters)

        new(basic, params_struct_array, order)
    end
end

DetailedTempoOutputResult() = DetailedTempoOutputResult(BasicTempoOutputResult(), FitParameter[])
DetailedTempoOutputResult(basic::BasicTempoOutputResult) = DetailedTempoOutputResult(basic, FitParameter[])

function Base.getproperty(dtr::DetailedTempoOutputResult, prop::Symbol)
    if prop == :basic
        return Base.getfield(dtr, :basic)
    elseif prop == :fit_parameters
        return Base.getfield(dtr, :fit_parameters)
    elseif prop == :fit_parameters_order
        return Base.getfield(dtr, :fit_parameters_order)
    elseif prop in fieldnames(BasicTempoOutputResult)
        return Base.getfield(dtr.basic, prop)
    elseif haskey(dtr.fit_parameters_order, prop)
        index = dtr.fit_parameters_order[prop]
        return dtr.fit_parameters[index]
    else
        error("Property $prop not found in DetailedTempoOutputResult.")
    end
end

# # Определим staticschema для DetailedTempoOutputResult
# StructArrays.staticschema(::Type{DetailedTempoOutputResult}) = NamedTuple{
#     (:basic, :fit_parameters, :fit_parameters_order, fieldnames(BasicTempoOutputResult)...),
#     Tuple{[fieldtype(DetailedTempoOutputResult, f) for f in fieldnames(DetailedTempoOutputResult)]..., 
#     [fieldtype(fieldtype(DetailedTempoOutputResult, :basic), a) for a in fieldnames(BasicTempoOutputResult)]...}
# }

# # Определим component для доступа к подполям basic и к полям fit_parameters
# function StructArrays.component(s::DetailedTempoOutputResult, key::Symbol)
#     if key in fieldnames(BasicTempoOutputResult)
#         return getproperty(getfield(s, :basic), key)
#     elseif haskey(s.fit_parameters_order, key)
#         index = s.fit_parameters_order[key]
#         return s.fit_parameters[index]
#     else
#         return getfield(s, key)
#     end
# end


function Base.show(io::IO, detailed::DetailedTempoOutputResult)
    indent = get(io, :indent, 0)
    println(io, ' '^indent, "Detailed Tempo Output Result:")
    show(IOContext(io, :indent => indent+4), detailed.basic)
    println(io, ' '^(indent+4), "Fit Parameters:")
    for param in detailed.fit_parameters
        show(IOContext(io, :indent => indent+8), param)
    end
end

struct TempoOutputError
    error_message::String
    error_details::String
    error_type::Symbol  # Например, :nan_values, :convergence_issue, :fit_error, etc.

    # Конструктор для пустой ошибки
    # TempoOutputError() = new("", "", :no_error)

    # Стандартный конструктор
    TempoOutputError(error_message::String, error_details::String, error_type::Symbol) = new(error_message, error_details, error_type)
end


# Конструктор с ключевыми словами и значениями по умолчанию
TempoOutputError(; 
    error_message::String = "",  
    error_details::String = "", 
    error_type::Symbol = :no_error
) = TempoOutputError(error_message, error_details, error_type)

function Base.show(io::IO, error::TempoOutputError)
    indent = get(io, :indent, 0)
    println(io, ' '^indent, "Tempo Output Error:")
    println(io, ' '^(indent+4), "Error message: ", error.error_message)
    println(io, ' '^(indent+4), "Error type: ", error.error_type)
end

# Результат одной итерации
struct InternalIterationTempoResult
    result::DetailedTempoOutputResult
    error::TempoOutputError
end

function Base.show(io::IO, int_iter::InternalIterationTempoResult)
    indent = get(io, :indent, 0)
    println(io, ' '^indent, "Internal iteration:")
    show(IOContext(io, :indent => indent+4), int_iter.result)
    show(IOContext(io, :indent => indent+4), int_iter.error)
end

InternalIterationTempoResult() = InternalIterationTempoResult(DetailedTempoOutputResult(), TempoOutputError())


# Общий результат, включая все итерации
struct GeneralTempoResult
    last_internal_iteration::InternalIterationTempoResult
    final_par_file::TempoParFile
    all_internal_iterations::Union{StructArray{InternalIterationTempoResult}, Nothing}
    all_global_iterations::Union{StructArray{GeneralTempoResult}, Nothing}

    # Конструктор для последней внутренней итерации
    GeneralTempoResult(last_internal_iteration::InternalIterationTempoResult, final_par_file::TempoParFile) = new(last_internal_iteration, final_par_file, nothing, nothing)

    # Конструктор для всех внутренних итераций
    function GeneralTempoResult(final_par_file::TempoParFile, all_internal_iterations::Vector{InternalIterationTempoResult})
        last_internal_iteration = isempty(all_internal_iterations) ? InternalIterationTempoResult() : all_internal_iterations[end]
        all_internal_iterations = isempty(all_internal_iterations) ? StructArray(InternalIterationTempoResult[]) : StructArray(all_internal_iterations)

        return new(last_internal_iteration, final_par_file, all_internal_iterations, nothing)
    end
    
    # Конструктор для всех глобальных итераций
    function GeneralTempoResult(all_global_iterations::Vector{GeneralTempoResult})
        last_internal_iteration = all_global_iterations[end].last_internal_iteration
        final_par_file = all_global_iterations[end].final_par_file
        all_internal_iterations = vcat([all_global_iterations[i].all_internal_iterations for i in 1:length(all_global_iterations)]...)

        return new(last_internal_iteration, final_par_file, all_internal_iterations, StructArray(all_global_iterations))
    end
end

function Base.show(io::IO, tempo_result::GeneralTempoResult)
    indent = get(io, :indent, 0)
    println(io, ' '^indent, "General Tempo Result:")
    show(IOContext(io, :indent => indent+4), tempo_result.last_internal_iteration)
    println(io, ' '^(indent+4), "Number of saved internal iterations: ", isnothing(tempo_result.all_internal_iterations) ? 0 : length(tempo_result.all_internal_iterations))
    println(io, ' '^(indent+4), "Number of saved global iterations: ",   isnothing(tempo_result.all_global_iterations)   ? 0 : length(tempo_result.all_global_iterations))
    # Здесь показываем поля структуры SingleTempoRunResult с соответствующими отступами
end

function extract_internal_iterations_values(all_iterations::StructArray{InternalIterationTempoResult}, param_name::Symbol, field_name::Symbol = :post_fit)
    # Проверка, существует ли поле в basic
    if param_name in fieldnames(BasicTempoOutputResult)
        return [getfield(itr.result.basic, param_name) for itr in all_iterations]
    end


    # Проверка, является ли параметр одним из фитируемых параметров
    if any(haskey(itr.result.fit_parameters_order, param_name) for itr in all_iterations)
        return [isempty(itr.result.fit_parameters)  ? NaN : getfield(itr.result.fit_parameters[itr.result.fit_parameters_order[param_name]], field_name) for itr in all_iterations]
    end

    # Если параметр не найден ни в basic, ни среди фитируемых параметров
    error("Parameter $param_name not found in basic or fit_parameters")
end

extract_internal_iterations_values(result::GeneralTempoResult, param_name::Symbol, field_name::Symbol = :post_fit) = extract_internal_iterations_values(result.all_internal_iterations, param_name, field_name)


#--------------------------------------------------------------------------------------------------------------

function parse_all_internal_interations_tempo_output(output::String, ::Type{Tempo2})::Vector{InternalIterationTempoResult}
    # Разделение на блоки итераций
    sections = split(output, "[tempo2.C:565] Complete fit")[2:end]  # Пропускаем первую секцию

    all_internal_interations = Vector{InternalIterationTempoResult}()

    for (niter, section) in enumerate(sections)
        section = String(section)
        internal_iteration_result = parse_internal_interation_tempo_output(section, Tempo2)

        # Проверка на наличие ошибок в секции

        # if internal_iteration_result.error !== TempoOutputError()
        #     break
        # end

        # Если ошибок нет, парсим результаты итерации

        push!(all_internal_interations, internal_iteration_result)
    end

    return all_internal_interations
end


function parse_internal_interation_tempo_output(section::String, ::Type{Tempo2})::InternalIterationTempoResult
    # Парсинг ошибок
    error = parse_tempo_output_error(section, Tempo2)

    # Получение BasicTempoOutputResult
    basic_result = parse_basic_tempo_output(section, Tempo2)

    if error !== TempoOutputError()
        # do nothing
    elseif basic_result.chisqr_red > 200.0
        error = TempoOutputError(error_type = :lost_connection)
        # return InternalIterationTempoResult(DetailedTempoOutputResult(basic_result), error)
    elseif basic_result.nfree < 1
        error = TempoOutputError(error_type = :not_enough_toas)
        # return InternalIterationTempoResult(DetailedTempoOutputResult(basic_result), error)
    # elseif basic_result.pre_post < 0.5 || basic_result.pre_post > 2.0
    #     error = TempoOutputError(error_type = :not_converged)
    #     # return InternalIterationTempoResult(DetailedTempoOutputResult(basic_result), error)
    end

    # Если обнаружена ошибка, возвращаем пустой DetailedTempoOutputResult и информацию об ошибке
    if error !== TempoOutputError()
        return InternalIterationTempoResult(DetailedTempoOutputResult(), error)
    end

    # Парсинг детальной информации о параметрах
    fit_parameters = parse_fit_parameters(section, Tempo2)

    # Создание DetailedTempoOutputResult
    detailed_result = DetailedTempoOutputResult(basic_result, fit_parameters)

    return InternalIterationTempoResult(detailed_result, error)
end

function parse_basic_tempo_output(section::String, ::Type{Tempo2})::BasicTempoOutputResult
    rms_regex        = r"RMS pre-fit residual = (\d+\.\d+|[-+]?[Nn][Aa][Nn]) \(us\), RMS post-fit residual = (\d+\.\d+|[-+]?[Nn][Aa][Nn]) \(us\)"
    rms_tn_regex     = r"RMS post-fit residual TN = (\d+\.\d+|[-+]?[Nn][Aa][Nn]) \(us\)"
    chisq_regex      = r"Fit Chisq = (\d+\.?\d*[eE]?[-+]?\d*|[-+]?[Nn][Aa][Nn])\s+Chisqr/nfree = (\d+\.?\d*[eE]?[-+]?\d*|[-+]?[Nn][Aa][Nn])/(0|\d+)\s*=\s*(\d+\.?\d*[eE]?[-+]?\d*|[-+]?[Nn][Aa][Nn])\s+pre/post = (\d+\.?\d*[eE]?[-+]?\d*|[-+]?[Nn][Aa][Nn])"
    fit_params_regex = r"Number of fit parameters: (\d+)"
    points_regex     = r"Number of points in fit = (\d+)"
    offset_regex     = r"Offset: ([\d\.\-eE]+|[-+]?[Nn][Aa][Nn]) ([\d\.\-eE]+|[-+]?[Nn][Aa][Nn]) offset_e\*sqrt\(n\) = ([\d\.\-eE]+|[-+]?[Nn][Aa][Nn]) n = (\d+)"

    # Функция для безопасного парсинга чисел, включая `-nan`
    function safe_parse_number(value::Union{Nothing, SubString{String}})::Float64
        if value === nothing || occursin(r"(?i)nan", value)
            return NaN
        else
            return parse(Float64, value)
        end
    end

    basic_result = BasicTempoOutputResult()

    # Найдем совпадения
    rms_match = match(rms_regex, section)
    rms_tn_match = match(rms_tn_regex, section)
    chisq_match = match(chisq_regex, section)
    fit_params_match = match(fit_params_regex, section)
    points_match = match(points_regex, section)
    offset_match = match(offset_regex, section)

    try
        basic_result = BasicTempoOutputResult(
            rms_match !== nothing ? safe_parse_number(rms_match[1]) : NaN,
            rms_match !== nothing ? safe_parse_number(rms_match[2]) : NaN,
            rms_tn_match !== nothing ? safe_parse_number(rms_tn_match[1]) : NaN,
            chisq_match !== nothing ? safe_parse_number(chisq_match[1]) : NaN,
            chisq_match !== nothing ? parse(Int, chisq_match[3]) : 0,
            chisq_match !== nothing ? safe_parse_number(chisq_match[4]) : NaN,
            chisq_match !== nothing ? safe_parse_number(chisq_match[5]) : NaN,
            fit_params_match !== nothing ? parse(Int, fit_params_match[1]) : 0,
            points_match !== nothing ? parse(Int, points_match[1]) : 0,
            offset_match !== nothing ? (safe_parse_number(offset_match[1]), safe_parse_number(offset_match[2])) : (NaN, NaN),
            offset_match !== nothing ? safe_parse_number(offset_match[3]) : NaN,
            NaN,
            NaN
        )
    catch error
        println("Ошибка при парсинге: $error")
    end

    return basic_result
end

function parse_fit_parameters(section::String, ::Type{Tempo2})::Vector{FitParameter}
    lines = split(section, '\n')
    param_block_indices = findall(x -> all(c -> c == '-', x) && x != "", lines)
    if length(param_block_indices) < 2
        return []  # Если нет двух дефисных линий, то возвращаем пустой список
    end
    param_lines = filter(x -> x[end] == 'Y' || x[end] == 'N', lines[param_block_indices[1]+1:param_block_indices[2]-1])

    regex = r"([\d\-\.eE]+)\s+([\d\-\.eE]+)\s+([\d\-\.eE]+)\s+([\d\-\.eE]+)\s+(Y|N)$"

    fit_parameters = []
    for line in param_lines
        parts = split(line)
        name = parts[1]
        match_result = match(regex, line)
        if match_result !== nothing
            pre_fit = parse(Float64, match_result[1])
            post_fit = parse(Float64, match_result[2])
            uncertainty = parse(Float64, match_result[3])
            difference = parse(Float64, match_result[4])
            fit_flag = match_result[5] == "Y"

            push!(fit_parameters, FitParameter(String(name), pre_fit, post_fit, uncertainty, difference, fit_flag))
        end
    end

    return fit_parameters
end


function parse_tempo_output_error(section::String, ::Type{Tempo2})::TempoOutputError
    # Примеры регулярных выражений для различных типов ошибок
    error_patterns = Dict(
        :nan_values => r"NaN",
        :convergence_issue => r"pre/post != 1",
        :fit_error => r"ERROR:.*"
        # Добавьте другие шаблоны по мере необходимости
    )

    for (error_type, pattern) in error_patterns
        match_result = match(pattern, section)
        if match_result !== nothing
            return TempoOutputError(
                error_message = String(match_result.match),
                error_details = "",  # Детали ошибки можно извлечь дополнительно
                error_type = error_type
            )
        end
    end

    return TempoOutputError()
end

#--------------------------------------------------------------------------------------------------------------



struct CalculatedResults
    # Поля для постфактум рассчитанных метрик, например:
    derived_metric::Float64
    # Другие рассчитанные метрики
end

struct TempoRunErrorOutput
    stderr_messages::String
    # Другие возможные поля для ошибок
end

struct SingleTempoRunResult{T <: AbstractTempoResult}
    tempo_output::T
    calculated_results::CalculatedResults
    error_output::TempoRunErrorOutput
    final_par_file::TempoParFile
    internal_iterations::Vector{T}
    global_iterations::Vector{SingleTempoRunResult{T}}

    # Конструктор для создания результата одной итерации (без подитераций)
    SingleTempoRunResult(
        tempo_output::T,
        calculated_results::CalculatedResults,
        error_output::TempoRunErrorOutput,
        final_par_file::TempoParFile
    ) where T <: AbstractTempoResult = new{T}(tempo_output, calculated_results, error_output, final_par_file, [], [])

    # Конструктор для создания результата с данными внутренних итераций                             
    SingleTempoRunResult(
        tempo_output::T,
        calculated_results::CalculatedResults,
        error_output::TempoRunErrorOutput,
        final_par_file::TempoParFile,
        internal_iterations::Vector{T}
    ) where T <: AbstractTempoResult = new{T}(tempo_output, calculated_results, error_output, final_par_file, internal_iterations, [])

    # Конструктор для создания результата с данными глобальных итераций
    SingleTempoRunResult(
        tempo_output::T,
        calculated_results::CalculatedResults,
        error_output::TempoRunErrorOutput,
        final_par_file::TempoParFile,
        global_iterations::Vector{SingleTempoRunResult{T}}
    ) where T <: AbstractTempoResult = new{T}(tempo_output, calculated_results, error_output, final_par_file, [], global_iterations)

    # Конструктор, который принимает массив результатов итераций
    SingleTempoRunResult(iteration_data::Vector{SingleTempoRunResult{T}}) where T <: AbstractTempoResult = new{T}(
        iteration_data[end].tempo_output,   # Используем результат последней итерации
        iteration_data[end].calculated_results,
        iteration_data[end].error_output,
        iteration_data[end].final_par_file,
        [],  # Пустой массив для внутренних итераций
        iteration_data  # Сохраняем весь массив итераций
    )
end

function Base.show(io::IO, run_result::SingleTempoRunResult)
    indent = get(io, :indent, 0)
    println(io, ' '^indent, "Single Tempo Run Result:")
    show(IOContext(io, :indent => indent+4), run_result.tempo_output)
    # Здесь показываем поля структуры SingleTempoRunResult с соответствующими отступами
end


# # Пример использования
# basic_result = BasicTempoOutputResult(...)  # Создание базового результата
# final_par_file_example = TempoParFile(...)  # Создание файла .par

# # Создание результата одной итерации
# single_result = SingleTempoRunResult(
#     basic_result,
#     CalculatedResults(...),  # Предполагаемые рассчитанные результаты
#     TempoRunErrorOutput("Some error occurred."),
#     final_par_file_example
# )

# # Создание результата с данными итераций
# iterations_results = [single_result, ...]  # Предполагаемый массив результатов итераций
# full_result = SingleTempoRunResult(
#     basic_result,
#     CalculatedResults(...),
#     TempoRunErrorOutput("Some error occurred."),
#     final_par_file_example,
#     iterations_results
# )

# # Предполагаем, что у нас есть массив результатов отдельных итераций с пустыми массивами итераций
# iteration_results = [SingleTempoRunResult(...), SingleTempoRunResult(...), ...]

# # Создаем общий результат, используя массив итераций
# final_result = SingleTempoRunResult(iteration_results)


