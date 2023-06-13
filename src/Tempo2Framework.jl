
function run_tempo2(par_file, tim_file, modified_tparams; add_flag="")
    chisqr = Inf
    modified_params_arr = []
    for (key,tparam) in modified_tparams
        push!(modified_params_arr, "-set", tparam.name, tparam.value)
    end
 
    try
        if add_flag == ""
            command = `tempo2 -f $par_file $tim_file $modified_params_arr`
        else
            command = `tempo2 -f $par_file $tim_file $(split(add_flag)) $modified_params_arr`
        end
#=
        if silent
            run(pipeline(command,stdout=devnull))
        else
            run(pipeline(command))
        end 
=#      
        out = readchomp(pipeline(command))
        out_split = split(out)
        i_chisqr = findall(x -> x == "Chisq", out_split)[end]
        chisqr_string = split(out_split[i_chisqr+5], "/")[1]
        chisqr = parse(Float64, chisqr_string)
    catch err
        println("Tempo2 failed")
        println(err)
        chisqr = Inf
    end
    return chisqr
end

function run_tempo2(par_file, tim_file; silent=true, add_flag="")
    chisqr = Inf

    try
        if add_flag == ""
            command = `tempo2 -f $par_file $tim_file`
        else
            command = `tempo2 -f $par_file $tim_file $(split(add_flag))`
        end   
        out = readchomp(pipeline(command))
        out_split = split(out)
        i_chisqr = findall(x -> x == "Chisq", out_split)[end]
        chisqr_string = split(out_split[i_chisqr+5], "/")[1]
        chisqr = parse(Float64, chisqr_string)
    catch err
        println("Tempo2 failed")
        println(err)
        chisqr = Inf
    end
    return chisqr
end

function calculate_t2!(tf::TempoFramework; add_refinement=0)

    par_file_init = tf.tsets.par_file_init
    par_file_work = "$(par_file_init[1:end-4])_work.par"
    par_file_out = "new.par"
    tim_file = tf.tsets.tim_file
    add_flag = tf.tsets.add_flag
    fit_XPBDOT = tf.tsets.fit_XPBDOT

    println("Generating derictories for parallel calculations: #workers = $(nworkers())")

    work_dir = pwd()
    for i in 2:nworkers()+1
        rm("./worker$i", force=true, recursive=true)
        mkdir("./worker$i")
        cp(tf.tsets.par_file_init, "./worker$i/" * tf.tsets.par_file_init, force=true)
        cp(tf.tsets.par_file_init, "./worker$i/" * par_file_work, force=true)
        cp(tf.tsets.tim_file, "./worker$i/" * tf.tsets.tim_file, force=true)
    end

    println("Obtaining a working parfile with fit_XPBDOT=$fit_XPBDOT")

    cp(par_file_init, par_file_work; force=true)
    if !haskey(tf.grid.params, :chisqr_min)
        tf.grid.params[:chisqr_min] = Inf
    end

    modified_tparams = Dict{String, TempoParameter}()
    if typeof(tf.test.alpha0) == Float64
        push!(modified_tparams, get_TempoParameter("alpha0", tf.test.alpha0))
    end
    if typeof(tf.test.beta0) == Float64
        push!(modified_tparams, get_TempoParameter("beta0", tf.test.beta0))
    end
    modified_tparams["eos"] = TempoParameter("EOS", tf.test.eosname)
    
    function get_ddstg_values_local(param1, param2; silent=true, only_keys = false, work_dir=work_dir)

        if myid() != 1
            cd(work_dir * "/worker$(myid())")
        end

        @printf "run %s = %10.6f, %s = %10.6f\n" tf.test.param1.name param1 tf.test.param2.name param2
  
        update_modifed_tparams!(modified_tparams, tf.tsets.params_first_step)
        modified_tparams["nits"] = TempoParameter("NITS", tf.tsets.nits_first_step)
        modified_tparams["gain"] = TempoParameter("GAIN", tf.tsets.gain_fisrt_step)
        push!(modified_tparams, get_TempoParameter(tf.test.param1.name, param1))
        push!(modified_tparams, get_TempoParameter(tf.test.param2.name, param2))
        modify_par_file(modified_tparams, par_file_work)

        chisqr = run_tempo2(par_file_work, tim_file; silent=silent, add_flag=add_flag)

        if tf.tsets.nits_second_step != 0
            update_modifed_tparams!(modified_tparams, tf.tsets.params_second_step)
            modified_tparams["nits"] = TempoParameter("NITS", tf.tsets.nits_second_step)
            modified_tparams["gain"] = TempoParameter("GAIN", tf.tsets.gain_second_step)
            modify_par_file(modified_tparams, par_file_out)
            chisqr = run_tempo2(par_file_out, tim_file; silent=silent, add_flag=add_flag)
        end

        tf.grid.params[:chisqr_min] = tf.grid.params[:chisqr_min] < chisqr ? tf.grid.params[:chisqr_min] : chisqr
#        temp_dict = read_params(Dict(:A1=>0.0, :E=>0.0, :T0=>0.0, :PB=>0.0, :OM=>0.0, :OMDOT=>0.0, :GAMMA=>0.0, :PBDOT=>0.0, :SINI=>0.0, :DTHETA=>0.0, :XDOT=>0.0, :DR=>0.0,:MA=>0.0, :MB =>0.0, :ALPHA0=>0.0, :BETA0=>0.0, :ALPHAA=>0.0, :BETAA=>0.0, :kA=>0.0), par_file_out)
    temp_dict = read_params(Dict(:A1=>0.0, :E=>0.0, :T0=>0.0, :PB=>0.0, :OM=>0.0, :OMDOT=>0.0, :GAMMA=>0.0, :PBDOT=>0.0, :SINI=>0.0, :H3 => 0.0, :VARSIGMA => 0.0, :DTHETA=>0.0, :XDOT=>0.0, :XPBDOT=>0.0, :DR=>0.0, :MTOT=>0.0, :M2 =>0.0, :ALPHA0=>0.0, :BETA0=>0.0, :ALPHAA=>0.0, :BETAA=>0.0, :kA=>0.0), par_file_out)

    @printf "DDSTG method m1 = %12.8f, m2 = %12.8f, χ2 = %8.3f, Δχ2 = %8.3f\n" temp_dict[:MTOT]-temp_dict[:M2] temp_dict[:M2] chisqr chisqr-tf.grid.params[:chisqr_min]

        ddstg_names = tuple(:chisqr, :eos_agn_chisqr, keys(temp_dict)...)
        ddstg_values = tuple(chisqr, chisqr, values(temp_dict)...)
        return (ddstg_names, ddstg_values)
    end

    chisqr_contours = tf.gsets.contours
    delta_chisqr_max = tf.gsets.delta_chisqr_max
    delta_chisqr_diff = tf.gsets.delta_chisqr_diff

    function niceplot_cell_selector(i_cell::Int64, j_cell::Int64, grid::Refinement2DGrid)
        chisqr_min = grid.params[:chisqr_min]
        chisqr_cell = @view grid.value[:chisqr][i_cell:i_cell+1,j_cell:j_cell+1]
        chisqr_cell_min = minimum(chisqr_cell)
        chisqr_cell_max = maximum(chisqr_cell)
        max_chisqr_case = (chisqr_cell_min < chisqr_min + delta_chisqr_max)
        diff_chisqr_case = (chisqr_cell_max - chisqr_cell_min > delta_chisqr_diff)
        contour_chisqr_case = any(chisqr_cell_min .< chisqr_min .+ chisqr_contours .< chisqr_cell_max)
        return max_chisqr_case && diff_chisqr_case || contour_chisqr_case
    end

    function contour_cell_selector(i_cell::Int64, j_cell::Int64, grid::Refinement2DGrid)
        chisqr_min = grid.params[:chisqr_min]
        chisqr_cell = @view grid.value[:chisqr][i_cell:i_cell+1,j_cell:j_cell+1]
        chisqr_cell_min = minimum(chisqr_cell)
        chisqr_cell_max = maximum(chisqr_cell)
#        min_chisqr_case = chisqr_cell_min <= chisqr_min < chisqr_cell_max
        contour_chisqr_case = any(chisqr_cell_min .< chisqr_min .+ chisqr_contours .< chisqr_cell_max)
        return contour_chisqr_case
    end

    function eos_agn_cell_selector(i_cell::Int64, j_cell::Int64, grid::Refinement2DGrid)
        chisqr_min = grid.params[:eos_agn_chisqr_min]
        chisqr_cell = @view grid.value[:eos_agn_chisqr][i_cell:i_cell+1,j_cell:j_cell+1]
        chisqr_cell_min = minimum(chisqr_cell)
        chisqr_cell_max = maximum(chisqr_cell)
        contour_chisqr_case = any(chisqr_cell_min .< chisqr_min .+ chisqr_contours .< chisqr_cell_max)
        return contour_chisqr_case && contour_cell_selector(i_cell, j_cell, grid)
    end

    if tf.gsets.refinement_type == "nice"
        sell_selector = niceplot_cell_selector
    elseif tf.gsets.refinement_type == "contour"
        sell_selector = contour_cell_selector
    elseif tf.gsets.refinement_type == "eos_agn"
        sell_selector = eos_agn_cell_selector
    end

    function calculate_params!(grid::Refinement2DGrid)
        if haskey(grid.params, :chisqr_min)
            grid.params[:chisqr_min] = min(minimum(grid.value[:chisqr]), grid.params[:chisqr_min])
        else
            grid.params[:chisqr_min] = minimum(grid.value[:chisqr])
        end
        if haskey(grid.params, :eos_agn_chisqr_min) && haskey(grid.value, :eos_agn_chisqr)
            grid.params[:eos_agn_chisqr_min] = min(minimum(grid.value[:eos_agn_chisqr]), grid.params[:eos_agn_chisqr_min])
        elseif haskey(grid.value, :eos_agn_chisqr)
            grid.params[:eos_agn_chisqr_min] = minimum(grid.value[:eos_agn_chisqr])
        end
        return nothing
    end

    if add_refinement == 0
        simple_parallel_precalculate_2DGrid(tf.grid, get_ddstg_values_local, calculate_params!)
        for i in 1:tf.gsets.N_refinement
            tf.grid = parallel_refine_2DGrid(tf.grid, get_ddstg_values_local, sell_selector, calculate_params!)
        end
    else
        if isempty(tf.grid.value)
            simple_parallel_precalculate_2DGrid(tf.grid, get_ddstg_values_local, calculate_params!)
        else
            for i in 1:add_refinement
                tf.grid = parallel_refine_2DGrid(tf.grid, get_ddstg_values_local, sell_selector, calculate_params!)
            end
            tf.gsets = GridSetttings(tf.gsets.N_refinement + add_refinement, tf.gsets.CL, tf.gsets.contours, tf.gsets.refinement_type, tf.gsets.delta_chisqr_max, tf.gsets.delta_chisqr_diff, tf.gsets.gr_in_chisqr)
        end
    end

    cut_ddstg_grid!(tf.grid)

    return tf
end