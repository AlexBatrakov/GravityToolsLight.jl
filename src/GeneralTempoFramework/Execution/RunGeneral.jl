#--------------------------------------------------------------------------------------------------------------
function run_tempo_single(bsets::BasicTempoSettings, gisets::Union{GlobalIterationsSettings,Nothing}=nothing)
    if gisets == nothing
        return run_tempo_basic(bsets)
    else
        return run_tempo_global_iters(bsets, gisets)
    end
end


function run_tempo_general(tf::GeneralTempoFramework; just_refine=false)
    work_dir = tf.tsets.basic_settings.work_dir
    cd(work_dir)
    
    if !just_refine
        rm("./runs", force=true, recursive=true)
        mkdir("./runs")
        GC.gc()  # ✅ Очистка после удаления и создания директорий
    end

    if tf.ref_sets.parallel
        for p in 1:nprocs()
            rm("./worker$p", force=true, recursive=true)
            mkdir("./worker$p")
            cp(tf.tsets.basic_settings.par_file_init, "$(work_dir)/worker$p/$(tf.tsets.basic_settings.par_file_init)", force=true)
            cp(tf.tsets.basic_settings.par_file_init, "$(work_dir)/worker$p/$(tf.tsets.basic_settings.par_file_init)", force=true)
            cp(tf.tsets.basic_settings.tim_file,      "$(work_dir)/worker$p/$(tf.tsets.basic_settings.tim_file)",      force=true)
        end
        GC.gc()  # ✅ Очистка после копирования файлов
    end

    function target_function(x, y, tf::GeneralTempoFramework=tf)
        work_dir = tf.tsets.basic_settings.work_dir

        if tf.ref_sets.parallel
            if myid() != 1
                work_dir = work_dir * "/worker$(myid())"
            end
        end

        x_name = tf.test_params.x.name
        y_name = tf.test_params.y.name

        bsets  = deepcopy(tf.tsets.basic_settings)
        gisets = deepcopy(tf.tsets.global_iter_settings)

        bsets.work_dir = work_dir

        if (tf.test_params.x.name != "TNRedAmp" && tf.test_params.x.name != "TNRedGam")
            push!(bsets.tparams, TP(tf.test_params.x.name, x, flag=-1))
        else
            push!(bsets.tparams, TP(tf.test_params.x.name, x))
        end
        if (tf.test_params.y.name != "TNRedAmp" && tf.test_params.y.name != "TNRedGam")
            push!(bsets.tparams, TP(tf.test_params.y.name, y, flag=-1))
        else
            push!(bsets.tparams, TP(tf.test_params.y.name, y))
        end
        println("Run started:  $x_name = $x, $y_name = $y")

        results = run_tempo_single(bsets, gisets)
        GC.gc()  # ✅ Очистка после выполнения `run_tempo_single`

        chisqr = results.last_internal_iteration.result.basic.chisqr
        rms_post_fit = results.last_internal_iteration.result.basic.rms_post_fit_residual_us
        pre_post = results.last_internal_iteration.result.basic.pre_post
        rms_tn_post_fit = results.last_internal_iteration.result.basic.rms_tn_post_fit_residual_us
        AD_objective = results.last_internal_iteration.result.basic.AD_objective
        chisqr_full = results.last_internal_iteration.result.basic.chisqr_full

        rm("$work_dir/design.matrix", force=true)
        cp(work_dir, "$(tf.tsets.basic_settings.work_dir)/runs/$x_name=$(x)_$y_name=$y", force=true)
        GC.gc()  # ✅ Очистка после копирования директорий

        values_arr = fill(NaN, length(tf.ref_sets.params_to_save))

        for (i, key) in enumerate(tf.ref_sets.params_to_save)
            if key == :chisqr
                values_arr[i] = results.last_internal_iteration.result.basic.chisqr
                values_arr[i] = isnan(values_arr[i]) ? Inf : values_arr[i]
            end
            if key == :rms_tn_post_fit
                values_arr[i] = results.last_internal_iteration.result.basic.rms_tn_post_fit_residual_us
            end
            if key == :rms_post_fit
                values_arr[i] = results.last_internal_iteration.result.basic.rms_post_fit_residual_us
            end
            if key == :pre_post
                values_arr[i] = results.last_internal_iteration.result.basic.pre_post
            end
            if key == :AD_objective
                values_arr[i] = results.last_internal_iteration.result.basic.AD_objective
            end
            if key == :chisqr_full
                values_arr[i] = results.last_internal_iteration.result.basic.chisqr_full
            end
            if haskey(results.final_par_file.tparams, key)
                values_arr[i] = Float64(results.final_par_file.tparams[key].value)
            end
        end

        GC.gc()  # ✅ Очистка после обработки параметров

        println("Run finished: $x_name = $x, $y_name = $y; chisqr = $chisqr, rms_post_fit = $rms_post_fit, rms_tn_post_fit = $rms_tn_post_fit, pre_post = $pre_post, AD_objective = $AD_objective, chisqr_full = $chisqr_full")

        return NamedTuple{tf.ref_sets.params_to_save}(values_arr)
    end

    function params_function!(grid::AdaptiveRefinement2DGrid)
    end

    if !just_refine
        tf.grid = calculate_2DGrid!(tf.grid, target_function, params_function!)
    else
        tf.grid = refine_2DGrid(tf.grid, target_function, params_function!)
    end

    GC.gc()  # ✅ Очистка после работы с сеткой

    return tf
end