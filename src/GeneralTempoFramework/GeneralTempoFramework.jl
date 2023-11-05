#--------------------------------------------------------------------------------------------------------------
# Определение абстрактного типа для версий Tempo
abstract type AbstractTempoVersion end

# Определение функции для всех типов Tempo
function get_tempo_command(::AbstractTempoVersion)
    error("This function must be implemented for each specific Tempo version.")
end

# Определение конструкторов и функций для Tempo
struct Tempo <: AbstractTempoVersion
    custom_data_directory::String
    Tempo() = new(get_tempo_directory()) # Конструктор по умолчанию
    Tempo(custom_dir::String) = new(custom_dir) # Конструктор с параметром
end

# Функция для получения пути к директории по умолчанию для Tempo
get_tempo_directory() = get(ENV, "TEMPO", "/default/path/to/tempo")

# Определение конструкторов и функций для Tempo2
struct Tempo2 <: AbstractTempoVersion
    custom_data_directory::String
    Tempo2() = new(get_tempo2_directory()) # Конструктор по умолчанию
    Tempo2(custom_dir::String) = new(custom_dir) # Конструктор с параметром
end

# Функция для получения пути к директории по умолчанию для Tempo2
get_tempo2_directory() = get(ENV, "TEMPO2", "/default/path/to/tempo2")

# Определение функций для получения команды
get_tempo_command(::Tempo) = "tempo"
get_tempo_command(::Tempo2) = "tempo2"


#--------------------------------------------------------------------------------------------------------------

struct TempoKeys
    silent::Bool
    print_output::Bool
    iterative_mode::Bool
    fit_EFACs_EQUADs::Bool
    TempoKeys(;silent = true, print_output = false, iterative_mode=true, fit_EFACs_EQUADs::Bool = false) = new(silent, print_output, iterative_mode, fit_EFACs_EQUADs)
end

function Base.show(io::IO, keys::TempoKeys)
    println(io, "Tempo keys:")
    print(io, "        Silent mode: ", keys.fit_EFACs_EQUADs)
    print(io, "        Fit EFAC and EQUAD parameters: ", keys.fit_EFACs_EQUADs)
	return nothing
end

#--------------------------------------------------------------------------------------------------------------
# Основные настройки Tempo
struct BasicTempoSettings{T <: AbstractTempoVersion}
    work_dir::String
    version::T
    par_file_init::String
    tim_file::String
    flags::String
    keys::TempoKeys
    tparams::Vector{GeneralTempoParameter}
end

function Base.show(io::IO, bsets::BasicTempoSettings)
    println(io, "Basic Tempo settings:")
    println(io, "   Working directory: ", bsets.work_dir)
    println(io, "   Version: ", bsets.version)
    println(io, "   Initial par file: ", bsets.par_file_init)
    println(io, "   Working tim file: ", bsets.tim_file)
    println(io, "   Selected additional flags: ", bsets.flags)
    println(io, "   ", bsets.keys)
    println(io, "   Tempo parameters used:", bsets.tparams)
	return nothing
end

BasicTempoSettings(;
    work_dir,
    version,
    par_file_init,
    tim_file,
    flags = "",
    keys = TempoKeys(),
    tparams = GeneralTempoParameter[]
    ) = BasicTempoSettings(
        work_dir,
        version,
        par_file_init,
        tim_file,
        flags,
        keys,
        tparams
        )

function run_tempo_basic(bsets::BasicTempoSettings)
    work_dir = bsets.work_dir  # Установка локальной переменной для директории
    cd(work_dir)
    
    # Создание объекта файла с начальными настройками
    par_file_init = TempoParFile("$(bsets.par_file_init)", new_name_suffix="_init")

    # Определение пути к файлу new.par и его удаление если он существует

    new_par_path = generate_par_file_path("new", "", work_dir)

    if isfile(new_par_path)
        rm(new_par_path)  # Удаление файла, чтобы избежать путаницы с предыдущими запусками
    end

    for tparam in bsets.tparams
        extend_par_file!(par_file_init, tparam)
    end

    write_par_file(par_file_init)

    command = `$(get_tempo_command(bsets.version)) -f $(par_file_init.name) $(bsets.tim_file) $([split(bsets.flags)...])`

        # Инициализация IOBuffer для захвата вывода
    output_io = IOBuffer()
    stderr_io = IOBuffer()

    # Запуск команды с перенаправлением вывода
    process = if bsets.keys.silent
        run(pipeline(command, stdout=output_io, stderr=stderr_io), wait=false)
    else
        run(pipeline(command, stdout=stdout, stderr=stderr_io), wait=false)
    end

    # Дождаться завершения процесса
    wait(process)

    # Считать вывод из IOBuffer
    output = String(take!(output_io))
    stderr_output = String(take!(stderr_io))

    # Сразу проверяем наличие ошибок в stderr_output и output
    if contains(stderr_output, "Ошибка") || contains(output, "NaN") || contains(output, "ERROR")
        println("Обнаружена ошибка в выводе. Ошибка выполнения $(typeof(bsets.version)).")
    end
    
    parsed_output = parse_tempo_output(output, bsets.version)

    # Если нужно, печатаем вывод в файл
    if bsets.keys.print_output
        output_filename = "$(par_file_init.name[1:end-4]).out"
        write(output_filename, output)
        # Если нужно, можно также сохранить stderr_output в файл
    end

    upd_par_path = generate_par_file_path(bsets.par_file_init, "upd", work_dir)


    # Проверяем существование файла
    if !isfile(new_par_path)
        println("Файл new_par_path не найден. Ошибка выполнения $(typeof(bsets.version)).")
    else
        cp(new_par_path, upd_par_path, force=true)
        par_file_upd = format_and_validate_par_file(upd_par_path, output, bsets)
    end

    return parsed_output, output, stderr_output
end


function format_and_validate_par_file(par_file_path::String, output::String, bsets::BasicTempoSettings)
    # Загрузка содержимого файла .par

    par_file = TempoParFile(par_file_path)

    # Проверка и добавление отсуствующих параметров
    # Ваша логика добавления параметров

    # тут добавляем отсуствующие глобальные параметры которые были при запуске. Если параметр присуствует то он либо отстался прежним либо был изменен
    for gparam in bsets.tparams
        if !haskey(par_file.tparams, gparam)
            extend_par_file!(par_file, gparam)
        end
    end

    # Форматирование файла
    # Ваша логика форматирования

    # Форматирование файла и добавление/обновление параметров EFAC и EQUAD, если это требуется
    if bsets.keys.fit_EFACs_EQUADs
        # Расчёт новых значений EFAC и EQUAD
        efacs_equads_params = calculate_EFACs_EQUADs(output, bsets)
        # Добавление/обновление EFAC и EQUAD в файле .par
        for efac_equad_param in efacs_equads_params
            extend_par_file!(par_file, efac_equad_param)
        end
    end

    # Сохранение отформатированного файла

    write_par_file(par_file)

    return par_file
end

# Функция для вызова парсера на основе типа версии
function parse_tempo_output(output, version::AbstractTempoVersion)
    parse_tempo_output(output, typeof(version))
end

# Определение функций для разных версий
function parse_tempo_output(output, ::Type{Tempo})
    # Реализация парсера для Tempo
end

function parse_tempo_output(output::String, ::Type{Tempo2})
    # Разделяем выходные данные на секции по итерациям
    sections = split(output, "Complete fit\n\n\n")
    
    # Убираем первый элемент, если он не содержит данных об итерации
    if !contains(sections[1], "RMS pre-fit residual")
        sections = sections[2:end]
    end


    # Создаем массив словарей для сохранения данных
    iterations_data = []
    
    for section in sections
        # Создаем словарь для сохранения данных текущей итерации
        results = Dict()
            
        # Извлечение информации о RMS pre-fit и post-fit residuals
        rms_regex = r"RMS pre-fit residual = (\d+\.\d+) \(us\), RMS post-fit residual = (\d+\.\d+) \(us\)"
        rms_match = match(rms_regex, section)
        if rms_match !== nothing
            results["RMS pre-fit residual (us)"] = parse(Float64, rms_match[1])
            results["RMS post-fit residual (us)"] = parse(Float64, rms_match[2])
        end
        
        # Обновляем регулярное выражение и процедуру извлечения информации о Fit Chisq, Chisqr/nfree и pre/post
        chisq_regex = r"Fit Chisq = (\d+\.?\d*[eE]?[-+]?\d*)\s+Chisqr/nfree = (\d+\.\d+)/(\d+) = (\d+\.\d+)\s+pre/post = (\d+(?:\.\d+)?)"
        chisq_match = match(chisq_regex, section)
        if chisq_match !== nothing
            results["Fit Chisq"] = parse(Float64, chisq_match[1])
            results["Chisqr"] = parse(Float64, chisq_match[2])
            results["nfree"] = parse(Int, chisq_match[3])
            results["Chisqr/nfree"] = parse(Float64, chisq_match[4])
            results["pre/post"] = parse(Float64, chisq_match[5])
        end
        
        # Извлечение информации о количестве параметров подгонки
        params_regex = r"Number of fit parameters: (\d+)"
        params_match = match(params_regex, section)
        if params_match !== nothing
            results["Number of fit parameters"] = parse(Int, params_match[1])
        end
        
        # Извлечение информации о количестве точек в подгонке
        points_regex = r"Number of points in fit = (\d+)"
        points_match = match(points_regex, section)
        if points_match !== nothing
            results["Number of points in fit"] = parse(Int, points_match[1])
        end
        
        # Извлечение информации об Offset
        offset_regex = r"Offset: ([\d\.\-eE]+) ([\d\.\-eE]+) offset_e\*sqrt\(n\) = ([\d\.\-eE]+) n = (\d+)"
        offset_match = match(offset_regex, section)
        if offset_match !== nothing
            results["Offset"] = (parse(Float64, offset_match[1]), parse(Float64, offset_match[2]))
            results["offset_e*sqrt(n)"] = parse(Float64, offset_match[3])
            results["n (number of points)"] = parse(Int, offset_match[4])
        end
        
        # Добавляем данные итерации в массив
        push!(iterations_data, results)

    end

    return iterations_data
end

#--------------------------------------------------------------------------------------------------------------        
# Настройки для глобальных итераций
struct GlobalIterationSettings
    iters::Int64
    nits::Vector{Int64}
    gain::Vector{Float64}
    tparams_local::Vector{Vector{GeneralTempoParameter}}
end

function Base.show(io::IO, gisets::GlobalIterationSettings)
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

GlobalIterationSettings(;
    iters,
    nits = 3 * ones(Int64, iters),
    gain = ones(Float64, iters),
    tparams_local = [Vector{GeneralTempoParameter}() for _ in 1:iters]
    ) = GlobalIterationSettings(
        iters,
        nits,
        gain,
        tparams_local
        )

function run_tempo_global_iter(bsets::BasicTempoSettings, gisets::GlobalIterationSettings)
    # Список для хранения результатов каждой итерации
    results = []

    # Удаление старых промежуточных .par файлов перед началом всех итераций
    for iter_file in readdir(bsets.work_dir)
        if occursin("_iter", iter_file) && (occursin(".par", iter_file) || occursin(".out", iter_file))
            rm(joinpath(bsets.work_dir, iter_file))  # Удаление файла
        end
    end

    # Используем `splitext` для получения имени файла без расширения
    base_par_name, ext = splitext(basename(bsets.par_file_init))

    # Если режим с итерациями, начните с начального .par файла
    previous_par_file = bsets.par_file_init

    # Определение пути к файлу new.par и его удаление если он существует
    new_par_path = joinpath(bsets.work_dir, "new.par")
    if isfile(new_par_path)
        rm(new_par_path)  # Удаление файла, чтобы избежать путаницы с предыдущими запусками
    end

    # Перебор каждой итерации
    for iter in 1:gisets.iters
        println("Итерация №$iter")
        
        # Обновление и применение настроек для текущей итерации
        current_flags = bsets.flags
        current_nits = iter <= length(gisets.nits) ? gisets.nits[iter] : gisets.nits[end]
        current_gain = iter <= length(gisets.gain) ? gisets.gain[iter] : gisets.gain[end]

        # Добавление NITS и GAIN в флаги, если они указаны
        current_flags *= current_nits > 0 ? " -set NITS $current_nits" : ""
        current_flags *= current_gain != 1.0 ? " -set GAIN $current_gain" : ""

        # Создание и обновление объекта настроек для итерации
    #    iter_bsets = deepcopy(bsets)

        current_par_file = joinpath(bsets.work_dir, base_par_name * "_iter$(iter)" * ext)

        if bsets.keys.iterative_mode
            # Если включен режим итераций, обновите iter_bsets для использования результатов предыдущей итерации
            cp(joinpath(bsets.work_dir, previous_par_file), current_par_file, force=true)
        else
            # В режиме независимых экспериментов всегда начинайте с исходного .par файла
            cp(joinpath(bsets.work_dir, bsets.par_file_init), current_par_file, force=true)
        end

#        iter_bsets.par_file_init = current_par_file

        iter_bsets = BasicTempoSettings(
        bsets.work_dir,
        bsets.version,
        current_par_file,
        bsets.tim_file,
        current_flags,
        bsets.keys,
        bsets.tparams
        )

#        iter_bsets.flags = current_flags

        # Обновляем локальные параметры для итерации, если они есть
        local_tparams = iter <= length(gisets.tparams_local) ? gisets.tparams_local[iter] : []

        # Добавляем или обновляем локальные параметры для текущей итерации
        for lparam in local_tparams
            found = false
            for (i, gparam) in enumerate(iter_bsets.tparams)
                if gparam.name == lparam.name
                    iter_bsets.tparams[i] = lparam  # Обновляем существующий параметр
                    found = true
                    break
                end
            end
            if !found
                push!(iter_bsets.tparams, lparam)  # Добавляем новый параметр
            end
        end

        # Запуск tempo с текущими настройками итерации
        parsed_output, output, stderr_output = run_tempo_basic(iter_bsets)

        # Проверка наличия ошибок и неправильных результатов TODO:stderr_output
        if bsets.keys.iterative_mode
            # Сразу проверяем наличие ошибок в stderr_output и output
            if contains(stderr_output, "Ошибка") || contains(output, "NaN") || contains(output, "ERROR")
                println("Обнаружена ошибка в выводе. Остановка итераций.")
                break
            end
        
            # Проверяем существование файла
            if !isfile(new_par_path)
                println("Файл new_par_path не найден. Остановка итераций.")
                break
            end
        end

        cp(new_par_path, current_par_file, force=true)

        # Форматирование и проверка файла .par после итерации

        formatted_par_file = format_and_validate_par_file(current_par_file, output, bsets)

        # Обновление ссылки на файл для следующей итерации и сохранение результатов
        if bsets.keys.iterative_mode
            # Только в режиме итераций обновляем previous_par_file
            previous_par_file = current_par_file
        end
        append!(results, parsed_output)
        println(parsed_output)
    end

    return results
end

#--------------------------------------------------------------------------------------------------------------
# Настройки для массива значений параметра (если у вас есть такой функционал)
struct ParameterSweepSettings
    param_name::String
    param_values::Vector{Float64}
    # Другие потенциальные поля...
end


#--------------------------------------------------------------------------------------------------------------
# Общие настройки Tempo, которые могут включать различные расширенные настройки
mutable struct GeneralTempoSettings{T <: AbstractTempoVersion}
    basic_settings::BasicTempoSettings{T}
    global_iter_settings::Union{GlobalIterationSettings, Nothing}
    param_sweep_settings::Union{ParameterSweepSettings, Nothing}
end

# Конструкторы
function GeneralTempoSettings(
    basic_settings::BasicTempoSettings{T}, 
    global_iter_settings::Union{GlobalIterationSettings, Nothing} = nothing,
    param_sweep_settings::Union{ParameterSweepSettings, Nothing} = nothing
    ) where {T <: AbstractTempoVersion}
    return GeneralTempoSettings{T}(basic_settings, global_iter_settings, param_sweep_settings)
end

#--------------------------------------------------------------------------------------------------------------






#--------------------------------------------------------------------------------------------------------------

mutable struct GeneralTempoFramework{T1 <: TestParameters, T2 <: RefinementSettings}
    tsets::GeneralTempoSettings
    test_params::T1
    ref_sets::T2
    grid::AdaptiveRefinement2DGrid
end

function Base.show(io::IO, tf::GeneralTempoFramework)
    println(io, "Tempo framework:")
    println(io, tf.tsets)
    println(io, tf.test_params)
    print(io, tf.ref_sets)
	return nothing
end


# function GeneralTempoFramework(test::GeneralTest, obs_params::ObsParams, gsets::GridSetttings)
#     param1_grid = collect(LinRange(test.param1.min, test.param1.max, test.param1.N))
#     param2_grid = collect(LinRange(test.param2.min, test.param2.max, test.param2.N))
#     grid = Refinement2DGrid(Dict(), param1_grid, param2_grid)
#     return PKFramework(test, obs_params, gsets, grid)
# end

GeneralTempoFramework(;tsets::GeneralTempoSettings, test_params::T1, ref_sets::T2) where {T1 <: TestParameters, T2 <: RefinementSettings} = TempoFramework{T1, T2}(tsets, test_params, ref_sets)

function GeneralTempoFramework(tsets::GeneralTempoSettings, test_params::TestParameters, ref_sets::RefinementSettings)
    grid = AdaptiveRefinement2DGrid(test_params.x, test_params.y, ref_sets)
    return GeneralTempoFramework(tsets, test_params, ref_sets, grid)
end


function calculate!(tf::GeneralTempoFramework)
    work_dir = tf.tsets.work_dir
    cd(work_dir)

    if tf.ref_sets.parallel
        for p in 1:nprocs()
            rm("./worker$p", force=true, recursive=true)
            mkdir("./worker$p")
            cp(tf.tsets.par_file_init, "$(work_dir)/worker$p/$(tf.tsets.par_file_init)", force=true)
            cp(tf.tsets.par_file_init, "$(work_dir)/worker$p/$(tf.tsets.par_file_init)", force=true)
            cp(tf.tsets.tim_file,      "$(work_dir)/worker$p/$(tf.tsets.tim_file)",      force=true)
        end
    end

    if tf.tsets.keys.fit_EFACs_EQUADs == true
        
    end


    function target_function(x, y, tf=tf)
        if tf.ref_sets.parallel
            p = myid()
            cd("$(work_dir)/worker$p")
        end

        par_file_init = TempoParFile("$(tf.tsets.par_file_init)", new_name_suffix="_init")
        
        x_name = tf.test_params.x.name_symbol
        y_name = tf.test_params.y.name_symbol

        par_file_init.tparams[x_name].value = x
        par_file_init.tparams[y_name].value = y

        write_par_file(par_file_init)

        for i in eachindex(tf.tsets.iters)

#        run_tempo_single(par_file_init, tim_file; silent=silent, add_flag=add_flag)
        end

        return (chisqr=sin(x^2*y)::Float64,)
    end

    function params_function!(grid::AdaptiveRefinement2DGrid)

    end

    calculate_2DGrid!(tf.grid, target_function, params_function!)

    return tf
end

function estimate_EFACs_EQUADs()
    return tparams
end
