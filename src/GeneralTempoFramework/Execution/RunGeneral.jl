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

    if tf.ref_sets.parallel
        for p in 1:nprocs()
            rm("./worker$p", force=true, recursive=true)
            mkdir("./worker$p")
            cp(tf.tsets.basic_settings.par_file_init, "$(work_dir)/worker$p/$(tf.tsets.basic_settings.par_file_init)", force=true)
            cp(tf.tsets.basic_settings.par_file_init, "$(work_dir)/worker$p/$(tf.tsets.basic_settings.par_file_init)", force=true)
            cp(tf.tsets.basic_settings.tim_file,      "$(work_dir)/worker$p/$(tf.tsets.basic_settings.tim_file)",      force=true)
        end
    end

    # if tf.tsets.keys.fit_EFACs_EQUADs == true
        
    # end

    function target_function(x, y, tf::GeneralTempoFramework=tf)

        lock_obj = ReentrantLock()

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

        push!(bsets.tparams, TP(tf.test_params.x.name, x, flag=0))
        push!(bsets.tparams, TP(tf.test_params.y.name, y, flag=0))

        lock(lock_obj) do
            println("Run started:  $x_name = $x, $y_name = $y")
        end

        results = run_tempo_single(bsets, gisets)

        chisqr = results.last_internal_iteration.result.basic.fit_chisq
        rms_post_fit = results.last_internal_iteration.result.basic.rms_post_fit_residual_us
        pre_post = results.last_internal_iteration.result.basic.pre_post
        rms_tn_post_fit = results.last_internal_iteration.result.basic.rms_tn_post_fit_residual_us

        lock(lock_obj) do
            println("Run finished: $x_name = $x, $y_name = $y; chisqr = $chisqr, rms_post_fit = $rms_post_fit, rms_tn_post_fit = $rms_tn_post_fit, pre_post = $pre_post")
        end

        # parfile_keys   = [:PB, :T0, :A1, :OM, :ECC, :PBDOT, :OMDOT, :M2, :MTOT, :GAMMA, :I, :IDOT, :AFAC, :BFAC, :AFACDOT, :BFACDOT, :AFAC2DOT, :BFAC2DOT, :P_LAMBDA, :P_ETA, :P_DELTA, :P_PHI]
        values_arr = Vector{Float64}(undef, length(tf.ref_sets.params_to_save))


        # temp_dict = Dict{Symbol, Float64}()

        # results.final_par_file.tparams[B].value

        for (i, key) in enumerate(tf.ref_sets.params_to_save)
            # parfile_values[i] = haskey(results.final_par_file.tparams, key) ? Float64(results.final_par_file.tparams[key].value) : NaN
            if key == :chisqr
                values_arr[i] = results.last_internal_iteration.result.basic.fit_chisq
            end
            if key == :rms_tn_post_fit
                values_arr[i] = results.last_internal_iteration.result.basic.rms_tn_post_fit_residual_us
            end
            if key == :rms_post_fit
                values_arr[i] = results.last_internal_iteration.result.basic.rms_post_fit_residual_us
            end
            if haskey(results.final_par_file.tparams, key)
                values_arr[i] = Float64(results.final_par_file.tparams[key].value)
            end
        end

        # return (chisqr=chisqr,)]
        # params_names  = tuple(:chisqr, keys(temp_dict)...)
        # params_values = tuple(chisqr, values(temp_dict)...)
        return NamedTuple{tf.ref_sets.params_to_save}(values_arr)
    end

    function params_function!(grid::AdaptiveRefinement2DGrid)

    end

    if !just_refine
        tf.grid = calculate_2DGrid!(tf.grid, target_function, params_function!)
    else
        tf.grid = refine_2DGrid(tf.grid, target_function, params_function!)
    end

    return tf
end
