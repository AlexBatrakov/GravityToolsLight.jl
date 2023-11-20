#--------------------------------------------------------------------------------------------------------------

abstract type AbstractTempoResult end

struct BasicTempoOutputResult <: AbstractTempoResult
    fit_chisq::Float64
    chisqr::Float64
    number_of_points_in_fit::Int
    number_of_fit_parameters::Int
    rms_pre_fit_residual_us::Float64
    nfree::Int
    offset::Tuple{Float64, Float64}
    offset_e_sqrt_n::Float64
    pre_post::Float64
    rms_post_fit_residual_us::Float64
    chisqr_nfree::Float64
end

function Base.show(io::IO, result::BasicTempoOutputResult)
    indent = get(io, :indent, 0)
    println(io, ' '^indent, "Basic Tempo Output Result:")
    println(io, ' '^(indent+4), "Fit Chisq: ", result.fit_chisq)
    println(io, ' '^(indent+4), "Chisqr: ", result.chisqr)
    println(io, ' '^(indent+4), "Number of Points in Fit: ", result.number_of_points_in_fit)
    println(io, ' '^(indent+4), "Number of Fit Parameters: ", result.number_of_fit_parameters)
    println(io, ' '^(indent+4), "RMS Pre-fit Residual (us): ", result.rms_pre_fit_residual_us)
    println(io, ' '^(indent+4), "Nfree: ", result.nfree)
    println(io, ' '^(indent+4), "Offset: (", result.offset[1], ", ", result.offset[2], ")")
    println(io, ' '^(indent+4), "Offset E * Sqrt(n): ", result.offset_e_sqrt_n)
    println(io, ' '^(indent+4), "Pre/Post: ", result.pre_post)
    println(io, ' '^(indent+4), "RMS Post-fit Residual (us): ", result.rms_post_fit_residual_us)
    println(io, ' '^(indent+4), "Chisqr/Nfree: ", result.chisqr_nfree)
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

# function Base.show(io::IO, param::FitParameter)
#     indent = get(io, :indent, 0)
#     println(io, ' '^indent, "Fit Parameter: ", param.name)
#     println(io, ' '^(indent+4), "Pre-fit: ", param.pre_fit)
#     println(io, ' '^(indent+4), "Post-fit: ", param.post_fit)
#     # и так далее для остальных полей
# end

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
    TempoOutputError() = new("", "", :no_error)

    # Стандартный конструктор
    TempoOutputError(error_message::String, error_details::String, error_type::Symbol) = new(error_message, error_details, error_type)
end

# Конструктор с ключевыми словами и значениями по умолчанию
TempoOutputError(; 
    error_message::String = "", 
    error_details::String = "", 
    error_type::Symbol = :unknown_error
) = TempoOutputError(iteration_number, error_message, error_details, error_type)


#--------------------------------------------------------------------------------------------------------------

function parse_all_interations_tempo_output(output::String, ::Type{Tempo2})::Tuple{Vector{DetailedTempoOutputResult}, Vector{TempoOutputError, Nothing}}
    # Разделение на блоки итераций
    sections = String.(split(output, "[tempo2.C:565] Complete fit")[2:end])  # Пропускаем первую секцию

    detailed_results = Vector{DetailedTempoOutputResult}()
    error_result = nothing

    for (niter, section) in enumerate(sections)

        detailed_result, error_in_section = parse_single_interation_tempo_output(section, niter, Tempo2)

        # Проверка на наличие ошибок в секции

        if error_in_section !== nothing
            error_in_section.iteration_number = niter
            error_result = error_in_section
            break
        end

        # Если ошибок нет, парсим результаты итерации

        push!(detailed_results, detailed_result)
    end

    return detailed_results, error_result
end


function parse_single_interation_tempo_output(section::String, niter::Int64, ::Type{Tempo2})::Tuple{DetailedTempoOutputResult, Union{TempoOutputError, Nothing}}
    # Парсинг ошибок
    error = parse_tempo_output_errors(section, niter, Tempo2)

    # Если обнаружена ошибка, возвращаем пустой DetailedTempoOutputResult и информацию об ошибке
    if error !== nothing
        return DetailedTempoOutputResult(), error
    end

    # Получение BasicTempoOutputResult
    basic_result = parse_basic_tempo_output(section, Tempo2)

    # Парсинг детальной информации о параметрах
    fit_parameters = parse_fit_parameters(section, Tempo2)

    # Создание DetailedTempoOutputResult
    detailed_result = DetailedTempoOutputResult(basic_result, fit_parameters)

    return detailed_result, nothing
end


function parse_basic_tempo_output(section::String, ::Type{Tempo2})::BasicTempoOutputResult

    # Создаем BasicTempoOutputResult для текущей итерации
    rms_regex = r"RMS pre-fit residual = (\d+\.\d+) \(us\), RMS post-fit residual = (\d+\.\d+) \(us\)"
    chisq_regex = r"Fit Chisq = (\d+\.?\d*[eE]?[-+]?\d*)\s+Chisqr/nfree = (\d+(?:\.\d*)?)(?:[eE][-+]?\d+)?/(\d+) = (\d+(?:\.\d*)?)(?:[eE][-+]?\d+)?\s+pre/post = (\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)"
    params_regex = r"Number of fit parameters: (\d+)"
    points_regex = r"Number of points in fit = (\d+)"
    offset_regex = r"Offset: ([\d\.\-eE]+) ([\d\.\-eE]+) offset_e\*sqrt\(n\) = ([\d\.\-eE]+) n = (\d+)"

    # Извлекаем значения
    rms_match = match(rms_regex, section)
    chisq_match = match(chisq_regex, section)
    params_match = match(params_regex, section)
    points_match = match(points_regex, section)
    offset_match = match(offset_regex, section)

    # Создаем объект BasicTempoOutputResult
    basic_result = BasicTempoOutputResult(
        parse(Float64, chisq_match[1]),  # fit_chisq
        parse(Float64, chisq_match[2]),  # chisqr
        parse(Int, points_match[1]),     # number_of_points_in_fit
        parse(Int, params_match[1]),     # number_of_fit_parameters
        parse(Float64, rms_match[1]),    # rms_pre_fit_residual_us
        parse(Int, chisq_match[3]),      # nfree
        (parse(Float64, offset_match[1]), parse(Float64, offset_match[2])),  # offset
        parse(Float64, offset_match[3]), # offset_e_sqrt_n
        parse(Float64, chisq_match[5]),  # pre_post
        parse(Float64, rms_match[2]),    # rms_post_fit_residual_us
        parse(Float64, chisq_match[4])   # chisqr_nfree
    )

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


function parse_tempo_output_errors(section::String, niter::Int64, ::Type{Tempo2})::Union{TempoOutputError, Nothing}
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
                iteration_number = niter,  # Необходимо определить номер итерации
                error_message = String(match_result.match),
                error_details = "",  # Детали ошибки можно извлечь дополнительно
                error_type = error_type
            )
        end
    end

    return nothing
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


