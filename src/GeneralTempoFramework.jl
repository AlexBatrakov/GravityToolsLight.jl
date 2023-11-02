@enum TempoVersion tempo=1 tempo2=2

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

mutable struct GeneralTempoSettings{T <: TempoVersion}
    work_dir::String
    version::T
    par_file_init::String
    tim_file::String
    flags::String
    keys::TempoKeys
    tparams::Vector{GeneralTempoParameter}
    iters::Int64
    nits::Vector{Int64}
    gain::Vector{Float64}
    tparams_local::Vector{Vector{GeneralTempoParameter}}
end

function Base.show(io::IO, tsets::GeneralTempoSettings)
    println(io, "Tempo settings:")
    println(io, "   Working directory: ", tsets.work_dir)
    println(io, "   Version: ", tsets.version)
    println(io, "   Initial par file: ", tsets.par_file_init)
    println(io, "   Working tim file: ", tsets.tim_file)
    println(io, "   Selected additional flags: ", tsets.flags)
    println(io, "   ", tsets.keys)
    println(io, "   Tempo parameters used:", tsets.tparams)
    println(io, "   Number of iterations: ", tsets.iters)
    for iter in 1:tsets.iters
        println(io, "   Step # $iter:")
        println(io, "       GAIN value: ", iter <= length(tsets.gain) ? tsets.gain[iter] : tsets.gain[end])
        println(io, "       NITS value: ", iter <= length(tsets.nits) ? tsets.nits[iter] : tsets.nits[end])
        println(io, "       Local tempo parameters: ", tsets.tparams_local[iter])
    end
	return nothing
end

GeneralTempoSettings(;
    work_dir,
    version,
    par_file_init,
    tim_file,
    flags = "",
    keys = TempoKeys(),
    tparams = GeneralTempoParameter[],
    iters,
    nits = 3 * ones(Int64, iters),
    gain = ones(Float64, iters),
    tparams_local = [Vector{GeneralTempoParameter}() for _ in 1:iters]
    ) = GeneralTempoSettings(
        work_dir,
        version,
        par_file_init,
        tim_file,
        flags,
        keys,
        tparams,
        iters,
        nits,
        gain,
        tparams_local
        )



function run_tempo_single(tsets::GeneralTempoSettings)
    work_dir = tsets.work_dir  # Установка локальной переменной для директории
    cd(work_dir)
    
    # Создание объекта файла с начальными настройками
    par_file_init = TempoParFile("$(tsets.par_file_init)", new_name_suffix="_init")

    # Определение пути к файлу new.par и его удаление если он существует

    new_par_path = generate_par_file_path("new", "", work_dir)

    if isfile(new_par_path)
        rm(new_par_path)  # Удаление файла, чтобы избежать путаницы с предыдущими запусками
    end

    for tparam in tsets.tparams
        extend_par_file!(par_file_init, tparam)
    end

    write_par_file(par_file_init)

    command = `$(tsets.version) -f $(par_file_init.name) $(tsets.tim_file) $([split(tsets.flags)...])`

        # Инициализация IOBuffer для захвата вывода
    output_io = IOBuffer()
    stderr_io = IOBuffer()

    # Запуск команды с перенаправлением вывода
    process = if tsets.keys.silent
        run(pipeline(command, stdout=output_io, stderr=stderr_io), wait=false)
    else
        run(pipeline(command, stdout=stdout, stderr=stderr_io), wait=false)
    end

    # Дождаться завершения процесса
    wait(process)

    # Считать вывод из IOBuffer
    output = String(take!(output_io))
    stderr_output = String(take!(stderr_io))

    parsed_results = parse_tempo_output(output)

    # Если нужно, печатаем вывод в файл
    if tsets.keys.print_output
        output_filename = "$(par_file_init.name[1:end-4]).out"
        write(output_filename, output)
        # Если нужно, можно также сохранить stderr_output в файл
    end

    upd_par_path = generate_par_file_path(tsets.par_file_init, "upd", work_dir)

    cp(new_par_path, upd_par_path, force=true)

    par_file_upd = format_and_validate_par_file(upd_par_path, output, tsets)

    return parsed_results, output, stderr_output
end


function format_and_validate_par_file(par_file_path::String, output::String, tsets::GeneralTempoSettings)
    # Загрузка содержимого файла .par

    par_file = TempoParFile(par_file_path)

    # Проверка и добавление отсуствующих параметров
    # Ваша логика добавления параметров

    # тут добавляем отсуствующие глобальные параметры которые были при запуске. Если параметр присуствует то он либо отстался прежним либо был изменен
    for gparam in tsets.tparams
        if !haskey(par_file.tparams, gparam)
            extend_par_file!(par_file, gparam)
        end
    end

    # Форматирование файла
    # Ваша логика форматирования

    # Форматирование файла и добавление/обновление параметров EFAC и EQUAD, если это требуется
    if tsets.keys.fit_EFACs_EQUADs
        # Расчёт новых значений EFAC и EQUAD
        efacs_equads_params = calculate_EFACs_EQUADs(output, tsets)
        # Добавление/обновление EFAC и EQUAD в файле .par
        for efac_equad_param in efacs_equads_params
            extend_par_file!(par_file, efac_equad_param)
        end
    end

    # Сохранение отформатированного файла

    write_par_file(par_file)

    return par_file
end

function parse_tempo_output(output::String)
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

        # Извлечение информации о Fit Chisq и Chisqr/nfree
        chisq_regex = r"Fit Chisq = (\d+\.?\d*[eE]?[-+]?\d*)\s+Chisqr/nfree = (\d+\.\d+)"
        chisq_match = match(chisq_regex, section)
        if chisq_match !== nothing
            results["Fit Chisq"] = parse(Float64, chisq_match[1])
            results["Chisqr/nfree"] = parse(Float64, chisq_match[2])
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

#using PyPlot

function plot_iterations_data(output_data)
    # Получаем ключи первой итерации для определения, какие параметры имеются
    params_keys = keys(output_data[1])
    
    # Определяем количество графиков
    num_plots = length(params_keys)
    
    # Создаем сетку графиков (пример для 4 графиков в ряду)
    num_cols = 4
    num_rows = cld(num_plots, num_cols)  # округление вверх для числа строк
    
    fig, axs = subplots(num_rows, num_cols, figsize=(15, num_rows*3))
    
    # Сделаем индексацию одномерной для удобства, если всего один график
    axs = reshape(axs, num_rows*num_cols)
    
    # Итерируем по параметрам и создаем графики
    for (index, param_key) in enumerate(params_keys)
        ax = axs[index]
        values = [iteration_data[param_key] for iteration_data in output_data]
        
        # Строим график
        ax.plot(values)
        ax.set_title(param_key)
        ax.set_xlabel("Iteration")
        ax.set_ylabel("Value")
        ax.grid(true)
        
        # Если достигли последнего ключа, прекращаем отрисовку
        if index == num_plots
            break
        end
    end
    
    # Если количество графиков не заполняет всю сетку, убираем пустые subplots
    for i in (num_plots+1):(num_rows*num_cols)
        fig.delaxes(axs[i])
    end
    
    # Показываем график
    tight_layout()
    show()
end

# Теперь вызываем функцию с данными
#plot_iterations_data(output_data)



function run_tempo_iterations(tsets::GeneralTempoSettings)
    # Список для хранения результатов каждой итерации
    results = []

    # Удаление старых промежуточных .par файлов перед началом всех итераций
    for iter_file in readdir(tsets.work_dir)
        if occursin("_iter", iter_file) && (occursin(".par", iter_file) || occursin(".out", iter_file))
            rm(joinpath(tsets.work_dir, iter_file))  # Удаление файла
        end
    end

    # Используем `splitext` для получения имени файла без расширения
    base_par_name, ext = splitext(basename(tsets.par_file_init))

    

    # Если режим с итерациями, начните с начального .par файла
    previous_par_file = tsets.par_file_init

    # Определение пути к файлу new.par и его удаление если он существует
    new_par_path = joinpath(tsets.work_dir, "new.par")
    if isfile(new_par_path)
        rm(new_par_path)  # Удаление файла, чтобы избежать путаницы с предыдущими запусками
    end

    # Перебор каждой итерации
    for iter in 1:tsets.iters
        println("Итерация №$iter")
        
        # Обновление и применение настроек для текущей итерации
        current_flags = tsets.flags
        current_nits = iter <= length(tsets.nits) ? tsets.nits[iter] : tsets.nits[end]
        current_gain = iter <= length(tsets.gain) ? tsets.gain[iter] : tsets.gain[end]

        # Добавление NITS и GAIN в флаги, если они указаны
        current_flags *= current_nits > 0 ? " -set NITS $current_nits" : ""
        current_flags *= current_gain != 1.0 ? " -set GAIN $current_gain" : ""

        # Создание и обновление объекта настроек для итерации
        iter_tsets = deepcopy(tsets)

        current_par_file = joinpath(tsets.work_dir, base_par_name * "_iter$(iter)" * ext)

        if tsets.keys.iterative_mode
            # Если включен режим итераций, обновите iter_tsets для использования результатов предыдущей итерации
            cp(joinpath(tsets.work_dir, previous_par_file), current_par_file, force=true)
        else
            # В режиме независимых экспериментов всегда начинайте с исходного .par файла
            cp(joinpath(tsets.work_dir, tsets.par_file_init), current_par_file, force=true)
        end

        iter_tsets.par_file_init = current_par_file


        iter_tsets.flags = current_flags

        # Обновляем локальные параметры для итерации, если они есть
        local_tparams = iter <= length(tsets.tparams_local) ? tsets.tparams_local[iter] : []

        # Добавляем или обновляем локальные параметры для текущей итерации
        for lparam in local_tparams
            found = false
            for (i, gparam) in enumerate(iter_tsets.tparams)
                if gparam.name == lparam.name
                    iter_tsets.tparams[i] = lparam  # Обновляем существующий параметр
                    found = true
                    break
                end
            end
            if !found
                push!(iter_tsets.tparams, lparam)  # Добавляем новый параметр
            end
        end

        # Запуск tempo с текущими настройками итерации
        output, stderr_output = run_tempo_single(iter_tsets)

        # Проверка наличия ошибок и неправильных результатов
        if tsets.keys.iterative_mode && (contains(stderr_output, "Ошибка") || !isfile(new_par_path) || contains(read(new_par_path, String), "NaN"))
            println("Обнаружена ошибка или неправильные результаты. Остановка итераций.")
            break
        end

        cp(new_par_path, current_par_file, force=true)

        # Форматирование и проверка файла .par после итерации

        formatted_par_file = format_and_validate_par_file(current_par_file, iter, tsets)

        # Обновление ссылки на файл для следующей итерации и сохранение результатов
        if tsets.keys.iterative_mode
            # Только в режиме итераций обновляем previous_par_file
            previous_par_file = current_par_file
            push!(results, (output, stderr_output))
        end
        push!(results, (output, stderr_output))
    end

    return results
end







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
