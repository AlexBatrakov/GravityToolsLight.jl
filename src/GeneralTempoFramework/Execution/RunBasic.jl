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
    
    # parsed_output = parse_tempo_output(output, bsets.version)
    # detailed_results, error_result = parse_all_internal_interations_tempo_output(output, typeof(bsets.version))
    all_internal_iterations = parse_all_internal_interations_tempo_output(output, typeof(bsets.version))

    # Если нужно, печатаем вывод в файл
    if bsets.keys.print_output
        output_filename = "$(par_file_init.name[1:end-4]).out"
        write(output_filename, output)
        # Если нужно, можно также сохранить stderr_output в файл
    end

    upd_par_path = generate_par_file_path(bsets.par_file_init, "upd", work_dir)
    par_file_upd = TempoParFile()

    # Проверяем существование файла
    if !isfile(new_par_path)
        println("Файл new_par_path не найден. Ошибка выполнения $(typeof(bsets.version)).")
    else
        cp(new_par_path, upd_par_path, force=true)
        par_file_upd = format_and_validate_par_file(upd_par_path, output, bsets)
    end

    general_tempo_result = GeneralTempoResult(par_file_upd, all_internal_iterations)
    return general_tempo_result

    # calculated_results = CalculatedResults(0.0)
    # error_output = TempoRunErrorOutput(stderr_output)
    # full_results = SingleTempoRunResult(
    #     detailed_results[end],
    #     calculated_results,
    #     error_output,
    #     par_file_upd,
    #     detailed_results
    # )
    # return full_results
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
        EFACs, EQUADs, log10EQUADs = calculate_EFACs_EQUADs(bsets)
        # Добавление/обновление EFAC и EQUAD в файле .par
        update_EFACs_EQUADs_in_par_file!(par_file, EFACs, log10EQUADs)
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

function parse_tempo_output(output::String, ::Type{Tempo2})::Vector{BasicTempoOutputResult}
    # Разделяем выходные данные на секции по итерациям
    sections = split(output, "Complete fit\n\n\n")
    
    # Убираем первый элемент, если он не содержит данных об итерации
    if !contains(sections[1], "RMS pre-fit residual")
        sections = sections[2:end]
    end

    # Создаем массив для сохранения данных
    iterations_data = BasicTempoOutputResult[]

    for section in sections
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
        result = BasicTempoOutputResult(
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

        push!(iterations_data, result)
    end

    return iterations_data
end


function parse_tempo_output_old(output::String, ::Type{Tempo2})
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
        chisq_regex = r"Fit Chisq = (\d+\.?\d*[eE]?[-+]?\d*)\s+Chisqr/nfree = (\d+(?:\.\d*)?)(?:[eE][-+]?\d+)?/(\d+) = (\d+(?:\.\d*)?)(?:[eE][-+]?\d+)?\s+pre/post = (\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)"
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