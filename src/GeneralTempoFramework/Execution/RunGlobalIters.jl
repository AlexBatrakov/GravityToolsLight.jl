

function run_tempo_global_iters(bsets::BasicTempoSettings, gisets::GlobalIterationsSettings)
    results_global = Vector{GeneralTempoResult}()

    # Удаление старых промежуточных .par файлов
    for iter_file in readdir(bsets.work_dir)
        if occursin("_iter", iter_file) && (occursin(".par", iter_file) || occursin(".out", iter_file))
            rm(joinpath(bsets.work_dir, iter_file))
        end
    end
    GC.gc()  # ✅ Очистка после удаления файлов

    base_par_name, ext = splitext(basename(bsets.par_file_init))
    previous_par_file = bsets.par_file_init
    new_par_path = joinpath(bsets.work_dir, "new.par")

    if isfile(new_par_path)
        rm(new_par_path)
    end
    GC.gc()  # ✅ Очистка после возможного удаления файла new.par

    for iter in 1:gisets.iters
        if gisets.print_output println("Global iteration №$iter started") end
        # println("While task. Used memory in GB: ", (Sys.total_memory() - Sys.free_memory()) / 1024^3)

        local_nits  = iter <= length(gisets.nits)  ? gisets.nits[iter]  : gisets.nits[end]
        local_gain  = iter <= length(gisets.gain)  ? gisets.gain[iter]  : gisets.gain[end]
        local_flags = iter <= length(gisets.flags) ? gisets.flags[iter] : gisets.flags[end]
        local_fit_EFACs_EQUADs = iter <= length(gisets.fit_EFACs_EQUADs) ? gisets.fit_EFACs_EQUADs[iter] : gisets.fit_EFACs_EQUADs[end]

        current_flags = bsets.flags
        current_flags *= local_nits > 0 ? " -set NITS $local_nits" : ""
        current_flags *= " -set GAIN $local_gain"
        current_flags *= local_flags != "" ? " $local_flags" : ""

        current_par_file = joinpath(bsets.work_dir, base_par_name * "_iter$(iter)" * ext)
        updated_current_par_file = joinpath(bsets.work_dir, base_par_name * "_iter$(iter)_upd" * ext)

        if gisets.keys.iterative_mode
            cp(joinpath(bsets.work_dir, previous_par_file), current_par_file, force=true)
        else
            cp(joinpath(bsets.work_dir, bsets.par_file_init), current_par_file, force=true)
        end

        iter_bsets = BasicTempoSettings(
            bsets.work_dir,
            bsets.version,
            current_par_file,
            bsets.tim_file,
            bsets.backends,
            current_flags,
            bsets.keys,
            bsets.tparams
        )
        bsets.keys.fit_EFACs_EQUADs = local_fit_EFACs_EQUADs

        local_tparams = iter <= length(gisets.tparams_local) ? gisets.tparams_local[iter] : []

        for lparam in local_tparams
            found = false
            for (i, gparam) in enumerate(iter_bsets.tparams)
                if gparam.name == lparam.name 
                    if iter == 1
                        iter_bsets.tparams[i].value = lparam.value != nothing ? lparam.value : iter_bsets.tparams[i].value
                        iter_bsets.tparams[i].flag  = lparam.flag
                    else
                        iter_bsets.tparams[i] = lparam
                    end
                    found = true
                    break
                end
            end
            if !found
                push!(iter_bsets.tparams, lparam)
            end
        end

        # ✅ Основной запуск TEMPO
        results_basic = run_tempo_basic(iter_bsets)

        # if isnan(results_basic.last_internal_iteration.result.basic.AD_objective) && iter > 1
        if local_fit_EFACs_EQUADs == false && iter > 1
            results_basic.last_internal_iteration.result.basic.AD_objective = results_global[end].last_internal_iteration.result.basic.AD_objective
            results_basic.last_internal_iteration.result.basic.chisqr_full  = results_global[end].last_internal_iteration.result.basic.chisqr_full
        end

        push!(results_global, results_basic)

        GC.gc()  # ✅ Очистка после run_tempo_basic

        results = results_basic.last_internal_iteration.result

        if gisets.print_output
            println("Global iteration №$iter finished\n" *
                    "   Local parameters: $(iter_bsets.tparams), fit EFACs & EQUADs: $(local_fit_EFACs_EQUADs)\n" *
                    "   chisqr = $(results.chisqr), nfree = $(results.nfree), chisqr/nfree = $(results.chisqr_red), RMS = $(results.rms_post_fit_residual_us), RMS_TN = $(results.rms_tn_post_fit_residual_us), Pre/post = $(results.pre_post), AD_objective = $(results.AD_objective), chisqr_full = $(results.chisqr_full)"
            )
        end

        if gisets.keys.iterative_mode
            if results_basic.last_internal_iteration.error !== TempoOutputError()
                println("Обнаружена ошибка. Остановка итераций.")
                GC.gc()  # ✅ Очистка перед выходом
                break
            end
        
            if !isfile(updated_current_par_file)
                println("Файл updated_current_par_file не найден. Остановка итераций.")
                GC.gc()  # ✅ Очистка перед выходом
                break
            end
        end

        try 
            cp(updated_current_par_file, current_par_file, force=true)
        catch error
            println("Файл updated_current_par_file не найден.")
        end

        if gisets.keys.iterative_mode
            previous_par_file = current_par_file
        end

        GC.gc()  # ✅ Очистка перед следующей итерацией
    end

    GC.gc()  # ✅ Очистка после всех итераций

    return GeneralTempoResult(results_global)
end