function run_tempo_basic(bsets::BasicTempoSettings)
    work_dir = bsets.work_dir
    cd(work_dir)
    
    par_file_init = TempoParFile("$(bsets.par_file_init)", new_name_suffix="_init")
    new_par_path = generate_par_file_path("new", "", work_dir)

    if isfile(new_par_path)
        rm(new_par_path)
    end

    for tparam in bsets.tparams
        extend_par_file!(par_file_init, tparam)
    end

    write_par_file(par_file_init)
    command = `$(get_tempo_command(bsets.version)) -f $(par_file_init.name) $(bsets.tim_file) $([split(bsets.flags)...])`

    # Инициализация IOBuffer
    output_io = IOBuffer()
    stderr_io = IOBuffer()

    # **Предварительная инициализация переменных для избежания UndefVarError**
    output = ""
    stderr_output = ""

    try
        # println("Запуск TEMPO2 в процессе $(myid()), worker $(getpid())")
        
        process = if bsets.keys.silent
            run(pipeline(command, stdout=output_io, stderr=stderr_io), wait=false)
        else
            run(pipeline(command, stdout=stdout, stderr=stderr_io), wait=false)
        end

        wait(process)

        output = String(take!(output_io))
        stderr_output = String(take!(stderr_io))

    finally
        close(output_io)
        close(stderr_io)
    end

    all_internal_iterations = parse_all_internal_interations_tempo_output(output, typeof(bsets.version))

    if isempty(all_internal_iterations)
        println("Ошибка выполнения $(typeof(bsets.version)): empty all_internal_iterations")
    elseif all_internal_iterations[end].error !== TempoOutputError()
        println("Ошибка выполнения $(typeof(bsets.version)): ", all_internal_iterations[end].error.error_type)
    end

    if bsets.keys.print_output
        output_filename = "$(par_file_init.name[1:end-4]).out"
        write(output_filename, output)
    end

    upd_par_path = generate_par_file_path(bsets.par_file_init, "upd", work_dir)
    par_file_upd = TempoParFile()

    if !isfile(new_par_path)
        println("Ошибка выполнения $(typeof(bsets.version)): файл $new_par_path не найден.")
    elseif all_internal_iterations[end].error == TempoOutputError()
        cp(new_par_path, upd_par_path, force=true)
        par_file_upd = format_and_validate_par_file(upd_par_path, output, bsets, all_internal_iterations)
    end

    general_tempo_result = GeneralTempoResult(par_file_upd, all_internal_iterations)

    GC.gc()

    return general_tempo_result
end


function format_and_validate_par_file(par_file_path::String, output::String, bsets::BasicTempoSettings, all_internal_iterations)
    # Загрузка содержимого файла .par

    par_file = TempoParFile(par_file_path)

    # Проверка и добавление отсуствующих параметров
    # Ваша логика добавления параметров

    # тут добавляем отсуствующие глобальные параметры которые были при запуске. Если параметр присуствует то он либо отстался прежним либо был изменен
    for gparam in bsets.tparams
        if !haskey(par_file.tparams, gparam.name_symbol)
            extend_par_file!(par_file, gparam)
        end
    end

    # Форматирование файла и добавление/обновление параметров EFAC и EQUAD, если это требуется
    if bsets.keys.fit_EFACs_EQUADs
        # Расчёт новых значений EFAC и EQUAD
        EFACs, EQUADs, log10EQUADs, AD_objective, chisqr_full = calculate_EFACs_EQUADs(bsets, time_start=par_file.tparams[:START].value, time_finish=par_file.tparams[:FINISH].value)
        # Добавление/обновление EFAC и EQUAD в файле .par
        all_internal_iterations[end].result.basic.AD_objective = AD_objective
        all_internal_iterations[end].result.basic.chisqr_full  = chisqr_full

        update_EFACs_EQUADs_in_par_file!(par_file, EFACs, log10EQUADs)
    end

    # Сохранение отформатированного файла

    write_par_file(par_file)

    return par_file
end