#--------------------------------------------------------------------------------------------------------------

function run_tempo_general!(tf::GeneralTempoFramework)
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
