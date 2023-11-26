

function run_tempo_global_iters(bsets::BasicTempoSettings, gisets::GlobalIterationsSettings)
    # Список для хранения результатов каждой итерации
    results_global = Vector{GeneralTempoResult}()

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
        println("Global iteration №$iter")
        
        # Обновление и применение настроек для текущей итерации
        current_flags = bsets.flags
        current_nits = iter <= length(gisets.nits) ? gisets.nits[iter] : gisets.nits[end]
        current_gain = iter <= length(gisets.gain) ? gisets.gain[iter] : gisets.gain[end]

        # Добавление NITS и GAIN в флаги, если они указаны
        current_flags *= current_nits > 0 ? " -set NITS $current_nits" : ""
        current_flags *= current_gain != 1.0 ? " -set GAIN $current_gain" : ""

        # Создание и обновление объекта настроек для итерации
        current_par_file = joinpath(bsets.work_dir, base_par_name * "_iter$(iter)" * ext)

        if gisets.keys.iterative_mode
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
        results_basic = run_tempo_basic(iter_bsets)
        push!(results_global, results_basic)


        # Проверка наличия ошибок и неправильных результатов TODO:stderr_output
        if gisets.keys.iterative_mode
            # Сразу проверяем наличие ошибок в stderr_output и output
            if results_basic.last_internal_iteration.error !== TempoOutputError()
                println("Обнаружена ошибка. Остановка итераций.")
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

        # formatted_par_file = format_and_validate_par_file(current_par_file, output, bsets)

        # Обновление ссылки на файл для следующей итерации и сохранение результатов
        if gisets.keys.iterative_mode
            # Только в режиме итераций обновляем previous_par_file
            previous_par_file = current_par_file
        end
        
    end

    if gisets.keys.save_global_iterations

    else
    
    end
    return GeneralTempoResult(results_global)
end