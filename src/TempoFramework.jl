
struct TempoSettings
    par_file_init::String
    tim_file::String
    add_flag::String
    fit_XPBDOT::Bool
    nits_first_step::Int64
    gain_fisrt_step::Float64
    params_first_step::Vector{Tuple{String, String, Int64}}
    nits_second_step::Int64
    gain_second_step::Float64
    params_second_step::Vector{Tuple{String, String, Int64}}
end

function Base.show(io::IO, tsets::TempoSettings)
    println(io, "Tempo settings:")
    println(io, "   Initial par file: ", tsets.par_file_init)
    println(io, "   Working tim file: ", tsets.tim_file)
    println(io, "   Selected additional flags: ", tsets.add_flag)
    println(io, "   Fit PBDOT to GR value: ", tsets.fit_XPBDOT)
    println(io, "   First step:")
    println(io, "   Number of iterations: ", tsets.nits_first_step)
    println(io, "   GAIN value for convergence stage: ", tsets.gain_fisrt_step)
    println(io, "   Parameters: ", tsets.params_first_step)
    println(io, "   Second step:")
    println(io, "   Number of iterations: ", tsets.nits_second_step)
    println(io, "   GAIN value for convergence stage: ", tsets.gain_second_step)
    print(io,   "   Parameters: ", tsets.params_second_step)
	return nothing
end

TempoSettings(;par_file_init, tim_file, add_flag, fit_XPBDOT, nits_first_step=5, gain_fisrt_step=1.0, params_first_step=Vector{Tuple{String, String, Int64}}[], nits_second_step=0, gain_second_step=1.0, params_second_step=Vector{Tuple{String, String, Int64}}[]) = TempoSettings(par_file_init, tim_file, add_flag, fit_XPBDOT, nits_first_step, gain_fisrt_step, params_first_step, nits_second_step, gain_second_step, params_second_step)


#struct STGTest <: AbstractGravityTest
#    psrname::String
#    eosname::String
#    log10alpha0::NamedTuple{(:min, :max, :N), Tuple{Float64, Float64, Int64}}
#    beta0::NamedTuple{(:min, :max, :N), Tuple{Float64, Float64, Int64}}
#end
#
#STGTest(;psrname, eosname, log10alpha0, beta0) = STGTest(psrname, eosname, log10alpha0, beta0)
#STGTest(;psrname, eosname, alpha0, beta0) = STGTest(psrname, eosname, (min = log10(abs(alpha0.min)), min = log10(abs(alpha0.max)), N = alpha0.N), beta0)
#
#function Base.show(io::IO, test::STGTest)
#    println(io, "STG test")
#    println(io, "Name of the pulsar: ", test.psrname)
#    println(io, "Equation of state: ", test.eosname)
#    println(io, "log10|alpha0| range: ", test.log10alpha0)
#    println(io, "beta0 range: ", test.beta0)
#	return nothing
#end
#
#struct MassMassTest <: AbstractGravityTest
#    psrname::String
#    eosname::String
#    alpha0::Float64
#    beta0::Float64
#    mpsr::NamedTuple{(:min, :max, :N), Tuple{Float64, Float64, Int64}}
#    mcomp::NamedTuple{(:min, :max, :N), Tuple{Float64, Float64, Int64}}
#    path_to_grids::String
#end
#
#MassMassTest(;psrname, eosname, alpha0, beta0, mpsr, mcomp, path_to_grids) = MassMassTest(psrname, eosname, alpha0, beta0, mpsr, mcomp, path_to_grids)
#
#function Base.show(io::IO, test::MassMassTest)
#    println(io, "Mass-mass test")
#    println(io, "Name of the pulsar: ", test.psrname)
#    println(io, "Equation of state: ", test.eosname)
#    println(io, "alpha0 value: ", test.alpha0)
#    println(io, "beta0 value: ", test.beta0)
#    println(io, "Pulsar mass range: ", test.mpsr)
#    println(io, "Companion mass range: ", test.mcomp)
#    println(io, "Path to grids: ", test.path_to_grids)
#	return nothing
#end

mutable struct TempoFramework{T <: AbstractGravityTest}
    test::T
    tsets::TempoSettings
    gsets::GridSetttings
    grid::SimpleGrid
end

function Base.show(io::IO, tf::TempoFramework)
    println(io, "Tempo framework:")
    println(io, tf.test)
    println(io, tf.tsets)
    print(io, tf.gsets)
	return nothing
end

#Base.copy(tf::TempoFramework) = TempoFramework(tf.test, tf.tsets, tf.gsets, tf.grid)

function TempoFramework(test::GeneralTest, obs_params::ObsParams, gsets::GridSetttings)
    param1_grid = collect(LinRange(test.param1.min, test.param1.max, test.param1.N))
    param2_grid = collect(LinRange(test.param2.min, test.param2.max, test.param2.N))
    grid = SimpleGrid(Dict(), param1_grid, param2_grid)
    return PKFramework(test, obs_params, gsets, grid)
end

TempoFramework(;test::T, tsets::TempoSettings, gsets::GridSetttings) where {T <: AbstractGravityTest} = TempoFramework{T}(test, tsets, gsets)

function TempoFramework(test::GeneralTest, tsets::TempoSettings, gsets::GridSetttings)
    param1_grid = collect(LinRange(test.param1.min, test.param1.max, test.param1.N))
    param2_grid = collect(LinRange(test.param2.min, test.param2.max, test.param2.N))
    grid = SimpleGrid(Dict(), param1_grid, param2_grid)
    return TempoFramework(test, tsets, gsets, grid)
end

#function TempoFramework(test::STGTest, tsets::TempoSettings, gsets::GridSetttings)
#    log10alpha0_grid = collect(LinRange(test.log10alpha0...))
#    beta0_grid = collect(LinRange(test.beta0...))
#    grid = SimpleGrid(Dict(), log10alpha0_grid, beta0_grid)
#    return TempoFramework(test, tsets, gsets, grid)
#end
#
#function TempoFramework(test::MassMassTest, tsets::TempoSettings, gsets::GridSetttings)
#    mpsr_grid = collect(LinRange(test.mpsr...))
#    mcomp_grid = collect(LinRange(test.mcomp...))
#    grid = SimpleGrid(Dict(), mpsr_grid, mcomp_grid)
#    return TempoFramework(test, tsets, gsets, grid)
#end

function write_new_params(params, par_file)

    open(par_file,"r") do file_par
        open("temp.par","w") do file_out
            for line in eachline(file_par)
                printed = false
                for param_name in keys(params)
                    if startswith(line, String(param_name))
                        if typeof(params[param_name]) <: Union{Float64,Int64}
                            println(file_out, param_name, "           ", @sprintf("%16s",params[param_name]))
                        elseif typeof(params[param_name]) == Bool
                            newline = replace(line, " 1 " => " $(params[param_name]*1) ", " 0 " => " $(params[param_name]*1) ")
                            println(file_out, newline)
                        elseif typeof(params[param_name]) == String
                            println(file_out, param_name, "           ", params[param_name])
                        end
                        printed = true
                    end
                end
                if printed == false
                    println(file_out, line)
                end
            end
        end
    end
    cp("temp.par", par_file; force=true)
end

function modify_par_file(params, par_file)
    function print_param(line, param_name, param_value)
        if typeof(param_value) <: Union{Float64,Int64}
            line = param_name * "           " * @sprintf("%16s",param_value)
        elseif typeof(param_value) == Bool
            line = replace(line, " 1 " => " $(param_value*1) ", " 0 " => " $(param_value*1) ")
        elseif typeof(param_value) == String
            line = param_name * "           " * param_value
        end
        return line
    end

    lines = readlines(par_file)
    for i in 1:length(lines)
        for param_name in keys(params)
            if startswith(lines[i], param_name)
                lines[i] = print_param(lines[i], param_name, params[param_name])
                delete!(params, param_name)
            end
        end
    end
    for param_name in keys(params)
        push!(lines, print_param("", param_name, params[param_name]))
    end
    open(par_file,"w") do file
        for line in lines
            println(file,line)
        end
    end
end

function read_chisqr()
    chisqr = Inf
    open("tempo.lis", "r") do f
        for line in eachline(f)
            if contains(line,r"Chisqr")
                try
                    chisqr = parse(Float64,line[15:24])
                    pre_post = parse(Float64,line[61:67])
                    if pre_post > 1.1
                        chisqr = Inf
                    end
                catch err
                    chisqr = Inf
                end
            end
        end
    end
    if isnan(chisqr)
        chisqr = Inf
    end
    return chisqr
end

function read_params(params, par_file)
    if isfile(par_file)
        open(par_file,"r") do file_out
            for line in eachline(par_file)
                for param_name in keys(params)
                    if startswith(line, String(param_name)*" ")
                        try
                            params[param_name] = parse(Float64,line[11:26])
#                        println(line)
                        catch err
                            params[param_name] = NaN
                        end
                    end
                end
            end
        end
    else
        for param_name in keys(params)
            params[param_name] = NaN
        end
    end
    return params
end

function run_tempo(par_file, tim_file; silent=true, add_flag="")
    chisqr = 1.0
    try
        if add_flag == ""
            command = `tempo -f $par_file $tim_file`
        else
            command = `tempo -f $par_file $tim_file $(split(add_flag))`
        end
        if silent
            run(pipeline(command,stdout=devnull))
        else
            run(pipeline(command))
        end
        chisqr = read_chisqr()
    catch err
        println("Tempo failed")
        println(err)
        chisqr = Inf
    end
    return chisqr
end

function read_chisqr_tempo_lis()
    lines = readlines("tempo.lis")
    chisqr = []
    for line in lines
        if contains(line,r"Chisqr")
            chisqr = append!(chisqr, parse(Float64,line[15:24]))
        end
    end
    return chisqr
end

function get_par_file_work(par_file_init, par_file_out, tim_file; add_flag="", fit_XPBDOT=true)
    par_file_work = "$(par_file_init[1:end-4])_work.par"
    cp(par_file_init, par_file_work; force=true)
    modify_par_file(Dict("ALPHA0"=>0.0, "BETA0"=>0.0, "NITS"=>30, "XPBDOT"=>fit_XPBDOT), par_file_work)
    run_tempo(par_file_work, tim_file; silent=true, add_flag="-c -no npulses.dat " * replace(add_flag, "-ni npulses.dat" => " "))
    cp(par_file_out, par_file_work; force=true)
    modify_par_file(Dict("NITS"=>5, "XPBDOT"=>false), par_file_work)
    chisqr = read_chisqr_tempo_lis()
    return mean(chisqr[10:30])
end

function get_tempo_format(name, value)
    if name == "alpha0"
        return "ALPHA0" => value
    elseif name == "log10alpha0"
        return "ALPHA0" => -exp10(value)
    elseif name == "beta0"
        return "BETA0" => value
    elseif name == "COSI"
        return "SINI" => sqrt(1-value^2)
    else
        return name => value
    end
end

#=
function calculate_old!(tf::TempoFramework; add_refinement=0)

    par_file_init = tf.tsets.par_file_init
    par_file_work = "$(par_file_init[1:end-4])_work.par"
    par_file_out = "$(tf.test.psrname).par"
    tim_file = tf.tsets.tim_file
    add_flag = tf.tsets.add_flag
    fit_XPBDOT = tf.tsets.fit_XPBDOT

    println("Obtaining a working parfile with fit_XPBDOT=$fit_XPBDOT")
    tf.grid.params[:chisqr_gr] = get_par_file_work(par_file_init, par_file_out, tim_file; add_flag=add_flag, fit_XPBDOT=fit_XPBDOT)
    tf.grid.params[:chisqr_min] = tf.gsets.gr_in_chisqr ? tf.grid.params[:chisqr_gr] : Inf

    modified_params = Dict()
    if typeof(tf.test.alpha0) == Float64
        push!(modified_params, get_tempo_format("alpha0", tf.test.alpha0))
    end
    if typeof(tf.test.beta0) == Float64
        push!(modified_params, get_tempo_format("beta0", tf.test.beta0))
    end
    
    function get_ddstg_values_local(param1, param2; silent=true)
        push!(modified_params, get_tempo_format(tf.test.param1.name, param1))
        push!(modified_params, get_tempo_format(tf.test.param2.name, param2))
        modified_params["EOS"] = tf.test.eosname
        modified_params["NITS"] = tf.tsets.nits
        modified_params["GAIN"] = 1.0
        modify_par_file(modified_params, par_file_work)
        chisqr = Inf
        if tf.tsets.all_conv == false
            chisqr = run_tempo(par_file_work, tim_file; silent=silent, add_flag=add_flag)
        end
        if isinf(chisqr) || tf.tsets.all_conv == true
            if tf.tsets.all_conv == false
                println("The first attempt failed. Reajustment of parameters.")
            end
            modified_params["NITS"] = tf.tsets.nits_conv
            modified_params["GAIN"] = tf.tsets.gain_conv
            modify_par_file(modified_params, par_file_work)
            chisqr = run_tempo(par_file_work, tim_file; silent=true, add_flag="-c "*add_flag)
            if !isinf(chisqr) && tf.tsets.nits != 0
                modified_params["NITS"] = tf.tsets.nits
                modified_params["GAIN"] = 1.0
                modify_par_file(modified_params, par_file_out)
                chisqr = run_tempo(par_file_out, tim_file; silent=true, add_flag=add_flag)
            end
        end
        tf.grid.params[:chisqr_min] = tf.grid.params[:chisqr_min] < chisqr ? tf.grid.params[:chisqr_min] : chisqr
        @printf "%s = %10.6f, %s = %10.6f, χ2 = %8.3f, Δχ2 = %8.3f\n" tf.test.param1.name param1 tf.test.param2.name param2 chisqr chisqr-tf.grid.params[:chisqr_min]
#        temp_dict = read_params(Dict(:A1=>0.0, :E=>0.0, :T0=>0.0, :PB=>0.0, :OM=>0.0, :OMDOT=>0.0, :GAMMA=>0.0, :PBDOT=>0.0, :SINI=>0.0, :DTHETA=>0.0, :XDOT=>0.0, :DR=>0.0,:MA=>0.0, :MB =>0.0, :ALPHA0=>0.0, :BETA0=>0.0, :ALPHAA=>0.0, :BETAA=>0.0, :kA=>0.0), par_file_out)
        temp_dict = read_params(Dict(:A1=>0.0, :E=>0.0, :T0=>0.0, :PB=>0.0, :OM=>0.0, :OMDOT=>0.0, :GAMMA=>0.0, :PBDOT=>0.0, :SINI=>0.0, :H3 => 0.0, :VARSIGMA => 0.0, :DTHETA=>0.0, :XDOT=>0.0, :XPBDOT=>0.0, :DR=>0.0, :MTOT=>0.0, :M2 =>0.0, :ALPHA0=>0.0, :BETA0=>0.0, :ALPHAA=>0.0, :BETAA=>0.0, :kA=>0.0), par_file_out)
        ddstg_names = tuple(:chisqr, keys(temp_dict)...)
        ddstg_values = tuple(chisqr, values(temp_dict)...)
        return (ddstg_names, ddstg_values)
    end

    chisqr_contours = tf.gsets.contours
    delta_chisqr_max = tf.gsets.delta_chisqr_max
    delta_chisqr_diff = tf.gsets.delta_chisqr_diff

    function niceplot_cell_selector(i_cell::Int64, j_cell::Int64, grid::SimpleGrid)
        chisqr_min = grid.params[:chisqr_min]
        chisqr_cell = @view grid.value[:chisqr][i_cell:i_cell+1,j_cell:j_cell+1]
        chisqr_cell_min = minimum(chisqr_cell)
        chisqr_cell_max = maximum(chisqr_cell)
        max_chisqr_case = (chisqr_cell_min < chisqr_min + delta_chisqr_max)
        diff_chisqr_case = (chisqr_cell_max - chisqr_cell_min > delta_chisqr_diff)
        contour_chisqr_case = any(chisqr_cell_min .< chisqr_min .+ chisqr_contours .< chisqr_cell_max)
        return max_chisqr_case && diff_chisqr_case || contour_chisqr_case
    end

    function contour_cell_selector(i_cell::Int64, j_cell::Int64, grid::SimpleGrid)
        chisqr_min = grid.params[:chisqr_min]
        chisqr_cell = @view grid.value[:chisqr][i_cell:i_cell+1,j_cell:j_cell+1]
        chisqr_cell_min = minimum(chisqr_cell)
        chisqr_cell_max = maximum(chisqr_cell)
#        min_chisqr_case = chisqr_cell_min <= chisqr_min < chisqr_cell_max
        contour_chisqr_case = any(chisqr_cell_min .< chisqr_min .+ chisqr_contours .< chisqr_cell_max)
        return contour_chisqr_case
    end

    if tf.gsets.refinement_type == "nice"
        sell_selector = niceplot_cell_selector
    elseif tf.gsets.refinement_type == "contour"
        sell_selector = contour_cell_selector
    end

    function calculate_params!(grid::SimpleGrid)
        if haskey(grid.params, :chisqr_min)
            grid.params[:chisqr_min] = min(minimum(grid.value[:chisqr]), grid.params[:chisqr_min])
        else
            grid.params[:chisqr_min] = minimum(grid.value[:chisqr])
        end
        return nothing
    end

    if add_refinement == 0
        precalculate_Grid(tf.grid, get_ddstg_values_local, calculate_params!)
        for i in 1:tf.gsets.N_refinement
            tf.grid = refine_Grid(tf.grid, get_ddstg_values_local, sell_selector, calculate_params!)
        end
    else
        for i in 1:add_refinement
            tf.grid = refine_Grid(tf.grid, get_ddstg_values_local, sell_selector, calculate_params!)
        end
        tf.gsets = GridSetttings(tf.gsets.N_refinement + add_refinement, tf.gsets.CL, tf.gsets.contours, tf.gsets.refinement_type, tf.gsets.delta_chisqr_max, tf.gsets.delta_chisqr_diff, tf.gsets.gr_in_chisqr)
    end

    cut_ddstg_grid!(tf.grid)
    return tf
end
=#

struct TempoParameter{T}
    name::String
    value::T
    flag::Int64
end

TempoParameter(name::String, value) = TempoParameter(name::String, value, -1)

function Base.show(io::IO, tparam::TempoParameter)
    print(io, "Tempo parameter: ")
    print(io, " ", tparam.name)
    print(io, " ", tparam.value)
    print(io, " ", tparam.flag)
	return nothing
end

function print_tparam(line::String, tparam::TempoParameter)
    flag_string = tparam.flag == -1 ? " " : "$(tparam.flag)"
    if typeof(tparam.value) <: Union{Float64,Int64}
        line = tparam.name * "           " * @sprintf("%16s",tparam.value) * " $(flag_string)"
    elseif tparam.value == "not changed"
        if length(split(line)) > 2
            line = replace(line, " 1 " => "  $(flag_string) ", " 0 " => " $(flag_string) ")
        else
            line = line * "  $(flag_string) "
        end
    elseif typeof(tparam.value) == String
        line = tparam.name * "           " * tparam.value
    end
    return line
end

function modify_par_file(tparams::Dict{String, TempoParameter}, par_file)
    lines = readlines(par_file)
    for i in 1:length(lines)
        for (name, tparam) in tparams
            if startswith(lines[i], tparam.name * " ")
                lines[i] = print_tparam(lines[i], tparam)
                delete!(tparams, name)
            end
        end
    end
    for tparam in values(tparams)
        push!(lines, print_tparam("", tparam))
    end
    open(par_file,"w") do file
        for line in lines
            println(file,line)
        end
    end
end

function get_TempoParameter(name, value, flag=-1)
    if name == "alpha0"
        return "alpha0" => TempoParameter("ALPHA0", value, flag)
    elseif name == "log10alpha0"
        return "alpha0" => TempoParameter("ALPHA0", -exp10(value), flag)
    elseif name == "beta0"
        return "beta0" => TempoParameter("BETA0", value, flag)
    elseif name == "mtot" || name == "m"
        return "mtot" => TempoParameter("MTOT", value, flag)
    elseif name == "mp" || name == "m2"
        return "m2" => TempoParameter("M2", value, flag)
    elseif name == "Pb"
        return "Pb" => TempoParameter("PB", value, flag)
    elseif name == "e0"
        return "e0" => TempoParameter("E", value, flag)
    elseif name == "x0"
        return "x0" => TempoParameter("A1", value, flag)
    elseif name == "omega0"
        return "omega0" => TempoParameter("OM", value, flag)
    elseif name == "COSI" || name == "cosi" 
        return "sini" => TempoParameter("SINI", sqrt(1 - value^2), flag)
    else
        return name => TempoParameter(name, value, flag)
    end
end

function update_pf_theory!(pf::DEFPhysicalFramework, name1, name2, param1, param2)
    key_change_theory = false
    if name1 == "log10alpha0"
        pf.theory.alpha0 = -exp10(param1)
        key_change_theory = true
    elseif  name1 == "alpha0"
        pf.theory.alpha0 = param1
        key_change_theory = true
    elseif name1 == "beta0"
        pf.theory.beta0  = param1
        key_change_theory = true
    end
    if name2 == "log10alpha0"
        pf.theory.alpha0 = -exp10(param2)
        key_change_theory = true
    elseif  name2 == "alpha0"
        pf.theory.alpha0 = param2
        key_change_theory = true
    elseif name2 == "beta0"
        pf.theory.beta0  = param2
        key_change_theory = true
    end
    if key_change_theory == true
        interpolate_mgrid!(pf)
    end
end

function update_modifed_tparams!(modified_tparams, params, m1, m2, Pb, e0, x0)
    for param in params
        if param[2] == "PK"
            if param[1] == "m2"
                push!(modified_tparams, get_TempoParameter("m2", m2, param[3]))
            elseif  param[1] == "mtot"
                push!(modified_tparams, get_TempoParameter("mtot", m1+m2, param[3]))
            elseif  param[1] == "Pb"
                push!(modified_tparams, get_TempoParameter("Pb", Pb, param[3]))
            elseif  param[1] == "e0"
                push!(modified_tparams, get_TempoParameter("e0", e0, param[3]))
            elseif  param[1] == "x0"
                push!(modified_tparams, get_TempoParameter("x0", x0, param[3]))
            end
        else
            push!(modified_tparams, get_TempoParameter(param...))
        end
    end
    return modified_tparams
end

function calculate!(tf::TempoFramework, pf::DEFPhysicalFramework, obs_params::ObsParams; add_refinement=0)

    par_file_init = tf.tsets.par_file_init
    par_file_work = "$(par_file_init[1:end-4])_work.par"
    par_file_out = "$(tf.test.psrname).par"
    tim_file = tf.tsets.tim_file
    add_flag = tf.tsets.add_flag
    fit_XPBDOT = tf.tsets.fit_XPBDOT

    pf.bnsys.K_params = obs_params.K
    pf.bnsys.psr.mass = obs_params.masses_init.m1
    pf.bnsys.comp.mass = obs_params.masses_init.m2
    pf.theory.alpha0  = 0.0
    pf.theory.beta0  = 0.0
    m1, m2 = find_best_masses(obs_params, pf)
    m1, m2, Pb, e0, x0 = optimize_PK_method(obs_params, pf)
    chisqr_PK_GR = GravityTools.get_chisqr(obs_params, pf)
    chisqr_PK_min = tf.gsets.gr_in_chisqr ? chisqr_PK_GR : Inf

    println("Obtaining a working parfile with fit_XPBDOT=$fit_XPBDOT")
    tf.grid.params[:chisqr_gr] = get_par_file_work(par_file_init, par_file_out, tim_file; add_flag=add_flag, fit_XPBDOT=fit_XPBDOT)
    if !haskey(tf.grid.params, :chisqr_min)
        tf.grid.params[:chisqr_min] = tf.gsets.gr_in_chisqr ? tf.grid.params[:chisqr_gr] : Inf
    end

    modified_tparams = Dict{String, TempoParameter}()
    if typeof(tf.test.alpha0) == Float64
        push!(modified_tparams, get_TempoParameter("alpha0", tf.test.alpha0))
        pf.theory.alpha0 = tf.test.alpha0
    end
    if typeof(tf.test.beta0) == Float64
        push!(modified_tparams, get_TempoParameter("beta0", tf.test.beta0))
        pf.theory.beta0 = tf.test.beta0
    end
    modified_tparams["eos"] = TempoParameter("EOS", tf.test.eosname)
    
    function get_ddstg_values_local(param1, param2; silent=true)
        @printf "run %s = %10.6f, %s = %10.6f\n" tf.test.param1.name param1 tf.test.param2.name param2

        pf.bnsys.K_params = obs_params.K
        pf.bnsys.psr.mass = obs_params.masses_init.m1
        pf.bnsys.comp.mass = obs_params.masses_init.m2
        update_pf_theory!(pf, tf.test.param1.name, tf.test.param2.name, param1, param2)
        m1, m2 = find_best_masses(obs_params, pf)
        m1, m2, Pb, e0, x0 = optimize_PK_method(obs_params, pf)
        chisqr_PK = get_chisqr(obs_params, pf)
        chisqr_PK_min = min(chisqr_PK_min, chisqr_PK)
        @printf "   PK method    m1 = %12.8f, m2 = %12.8f, χ2 = %8.3f, Δχ2 = %8.3f\n" m1 m2 chisqr_PK chisqr_PK-chisqr_PK_min
    
        update_modifed_tparams!(modified_tparams, tf.tsets.params_first_step, m1, m2, Pb, e0, x0)
        modified_tparams["nits"] = TempoParameter("NITS", tf.tsets.nits_first_step)
        modified_tparams["gain"] = TempoParameter("GAIN", tf.tsets.gain_fisrt_step)
        push!(modified_tparams, get_TempoParameter(tf.test.param1.name, param1))
        push!(modified_tparams, get_TempoParameter(tf.test.param2.name, param2))
        modify_par_file(modified_tparams, par_file_work)

        chisqr = run_tempo(par_file_work, tim_file; silent=silent, add_flag=add_flag)

        if tf.tsets.nits_second_step != 0
            update_modifed_tparams!(modified_tparams, tf.tsets.params_second_step, m1, m2, Pb, e0, x0)
            modified_tparams["nits"] = TempoParameter("NITS", tf.tsets.nits_second_step)
            modified_tparams["gain"] = TempoParameter("GAIN", tf.tsets.gain_second_step)
            modify_par_file(modified_tparams, par_file_out)
            chisqr = run_tempo(par_file_out, tim_file; silent=silent, add_flag=add_flag)
        end

        tf.grid.params[:chisqr_min] = tf.grid.params[:chisqr_min] < chisqr ? tf.grid.params[:chisqr_min] : chisqr
#        temp_dict = read_params(Dict(:A1=>0.0, :E=>0.0, :T0=>0.0, :PB=>0.0, :OM=>0.0, :OMDOT=>0.0, :GAMMA=>0.0, :PBDOT=>0.0, :SINI=>0.0, :DTHETA=>0.0, :XDOT=>0.0, :DR=>0.0,:MA=>0.0, :MB =>0.0, :ALPHA0=>0.0, :BETA0=>0.0, :ALPHAA=>0.0, :BETAA=>0.0, :kA=>0.0), par_file_out)
    temp_dict = read_params(Dict(:A1=>0.0, :E=>0.0, :T0=>0.0, :PB=>0.0, :OM=>0.0, :OMDOT=>0.0, :GAMMA=>0.0, :PBDOT=>0.0, :SINI=>0.0, :H3 => 0.0, :VARSIGMA => 0.0, :DTHETA=>0.0, :XDOT=>0.0, :XPBDOT=>0.0, :DR=>0.0, :MTOT=>0.0, :M2 =>0.0, :ALPHA0=>0.0, :BETA0=>0.0, :ALPHAA=>0.0, :BETAA=>0.0, :kA=>0.0), par_file_out)

    @printf "   DDSTG method m1 = %12.8f, m2 = %12.8f, χ2 = %8.3f, Δχ2 = %8.3f\n" temp_dict[:MTOT]-temp_dict[:M2] temp_dict[:M2] chisqr chisqr-tf.grid.params[:chisqr_min]

        ddstg_names = tuple(:chisqr, keys(temp_dict)...)
        ddstg_values = tuple(chisqr, values(temp_dict)...)
        return (ddstg_names, ddstg_values)
    end

    chisqr_contours = tf.gsets.contours
    delta_chisqr_max = tf.gsets.delta_chisqr_max
    delta_chisqr_diff = tf.gsets.delta_chisqr_diff

    function niceplot_cell_selector(i_cell::Int64, j_cell::Int64, grid::SimpleGrid)
        chisqr_min = grid.params[:chisqr_min]
        chisqr_cell = @view grid.value[:chisqr][i_cell:i_cell+1,j_cell:j_cell+1]
        chisqr_cell_min = minimum(chisqr_cell)
        chisqr_cell_max = maximum(chisqr_cell)
        max_chisqr_case = (chisqr_cell_min < chisqr_min + delta_chisqr_max)
        diff_chisqr_case = (chisqr_cell_max - chisqr_cell_min > delta_chisqr_diff)
        contour_chisqr_case = any(chisqr_cell_min .< chisqr_min .+ chisqr_contours .< chisqr_cell_max)
        return max_chisqr_case && diff_chisqr_case || contour_chisqr_case
    end

    function contour_cell_selector(i_cell::Int64, j_cell::Int64, grid::SimpleGrid)
        chisqr_min = grid.params[:chisqr_min]
        chisqr_cell = @view grid.value[:chisqr][i_cell:i_cell+1,j_cell:j_cell+1]
        chisqr_cell_min = minimum(chisqr_cell)
        chisqr_cell_max = maximum(chisqr_cell)
#        min_chisqr_case = chisqr_cell_min <= chisqr_min < chisqr_cell_max
        contour_chisqr_case = any(chisqr_cell_min .< chisqr_min .+ chisqr_contours .< chisqr_cell_max)
        return contour_chisqr_case
    end

    if tf.gsets.refinement_type == "nice"
        sell_selector = niceplot_cell_selector
    elseif tf.gsets.refinement_type == "contour"
        sell_selector = contour_cell_selector
    end

    function calculate_params!(grid::SimpleGrid)
        if haskey(grid.params, :chisqr_min)
            grid.params[:chisqr_min] = min(minimum(grid.value[:chisqr]), grid.params[:chisqr_min])
        else
            grid.params[:chisqr_min] = minimum(grid.value[:chisqr])
        end
        return nothing
    end

    if add_refinement == 0
        precalculate_Grid(tf.grid, get_ddstg_values_local, calculate_params!)
        for i in 1:tf.gsets.N_refinement
            tf.grid = refine_Grid(tf.grid, get_ddstg_values_local, sell_selector, calculate_params!)
        end
    else
        for i in 1:add_refinement
            tf.grid = refine_Grid(tf.grid, get_ddstg_values_local, sell_selector, calculate_params!)
        end
        tf.gsets = GridSetttings(tf.gsets.N_refinement + add_refinement, tf.gsets.CL, tf.gsets.contours, tf.gsets.refinement_type, tf.gsets.delta_chisqr_max, tf.gsets.delta_chisqr_diff, tf.gsets.gr_in_chisqr)
    end

    cut_ddstg_grid!(tf.grid)

    return tf, pf, obs_params
end

function cut_ddstg_grid!(grid::SimpleGrid)
    chisqr_min = grid.params[:chisqr_min]
    chisqr_max_mesuared = maximum(filter(x->x<Inf, grid.value[:chisqr]))
    grid.value[:chisqr_cut] = deepcopy(grid.value[:chisqr])
    for i in 1:grid.N_x, j in 1:grid.N_y
        if grid.value[:chisqr][i,j] > chisqr_max_mesuared
            grid.value[:chisqr_cut][i,j] = chisqr_max_mesuared
        end
    end
end

function update_modifed_tparams!(modified_tparams, params)
    for param in params
        push!(modified_tparams, get_TempoParameter(param...))
    end
    return modified_tparams
end

function calculate!(tf::TempoFramework; add_refinement=0)

    par_file_init = tf.tsets.par_file_init
    par_file_work = "$(par_file_init[1:end-4])_work.par"
    par_file_out = "$(tf.test.psrname).par"
    tim_file = tf.tsets.tim_file
    add_flag = tf.tsets.add_flag
    fit_XPBDOT = tf.tsets.fit_XPBDOT

    println("Obtaining a working parfile with fit_XPBDOT=$fit_XPBDOT")


    modified_tparams = Dict{String, TempoParameter}()
    if typeof(tf.test.alpha0) == Float64
        push!(modified_tparams, get_TempoParameter("alpha0", tf.test.alpha0))
    end
    if typeof(tf.test.beta0) == Float64
        push!(modified_tparams, get_TempoParameter("beta0", tf.test.beta0))
    end
    modified_tparams["eos"] = TempoParameter("EOS", tf.test.eosname)
    
    function get_ddstg_values_local(param1, param2; silent=true)
        @printf "run %s = %10.6f, %s = %10.6f\n" tf.test.param1.name param1 tf.test.param2.name param2
  
        update_modifed_tparams!(modified_tparams, tf.tsets.params_first_step)
        modified_tparams["nits"] = TempoParameter("NITS", tf.tsets.nits_first_step)
        modified_tparams["gain"] = TempoParameter("GAIN", tf.tsets.gain_fisrt_step)
        push!(modified_tparams, get_TempoParameter(tf.test.param1.name, param1))
        push!(modified_tparams, get_TempoParameter(tf.test.param2.name, param2))
        modify_par_file(modified_tparams, par_file_work)

        chisqr = run_tempo(par_file_work, tim_file; silent=silent, add_flag=add_flag)

        if tf.tsets.nits_second_step != 0
            update_modifed_tparams!(modified_tparams, tf.tsets.params_second_step)
            modified_tparams["nits"] = TempoParameter("NITS", tf.tsets.nits_second_step)
            modified_tparams["gain"] = TempoParameter("GAIN", tf.tsets.gain_second_step)
            modify_par_file(modified_tparams, par_file_out)
            chisqr = run_tempo(par_file_out, tim_file; silent=silent, add_flag=add_flag)
        end

        tf.grid.params[:chisqr_min] = tf.grid.params[:chisqr_min] < chisqr ? tf.grid.params[:chisqr_min] : chisqr
#        temp_dict = read_params(Dict(:A1=>0.0, :E=>0.0, :T0=>0.0, :PB=>0.0, :OM=>0.0, :OMDOT=>0.0, :GAMMA=>0.0, :PBDOT=>0.0, :SINI=>0.0, :DTHETA=>0.0, :XDOT=>0.0, :DR=>0.0,:MA=>0.0, :MB =>0.0, :ALPHA0=>0.0, :BETA0=>0.0, :ALPHAA=>0.0, :BETAA=>0.0, :kA=>0.0), par_file_out)
    temp_dict = read_params(Dict(:A1=>0.0, :E=>0.0, :T0=>0.0, :PB=>0.0, :OM=>0.0, :OMDOT=>0.0, :GAMMA=>0.0, :PBDOT=>0.0, :SINI=>0.0, :H3 => 0.0, :VARSIGMA => 0.0, :DTHETA=>0.0, :XDOT=>0.0, :XPBDOT=>0.0, :DR=>0.0, :MTOT=>0.0, :M2 =>0.0, :ALPHA0=>0.0, :BETA0=>0.0, :ALPHAA=>0.0, :BETAA=>0.0, :kA=>0.0), par_file_out)

    @printf "   DDSTG method m1 = %12.8f, m2 = %12.8f, χ2 = %8.3f, Δχ2 = %8.3f\n" temp_dict[:MTOT]-temp_dict[:M2] temp_dict[:M2] chisqr chisqr-tf.grid.params[:chisqr_min]

        ddstg_names = tuple(:chisqr, keys(temp_dict)...)
        ddstg_values = tuple(chisqr, values(temp_dict)...)
        return (ddstg_names, ddstg_values)
    end

    chisqr_contours = tf.gsets.contours
    delta_chisqr_max = tf.gsets.delta_chisqr_max
    delta_chisqr_diff = tf.gsets.delta_chisqr_diff

    function niceplot_cell_selector(i_cell::Int64, j_cell::Int64, grid::SimpleGrid)
        chisqr_min = grid.params[:chisqr_min]
        chisqr_cell = @view grid.value[:chisqr][i_cell:i_cell+1,j_cell:j_cell+1]
        chisqr_cell_min = minimum(chisqr_cell)
        chisqr_cell_max = maximum(chisqr_cell)
        max_chisqr_case = (chisqr_cell_min < chisqr_min + delta_chisqr_max)
        diff_chisqr_case = (chisqr_cell_max - chisqr_cell_min > delta_chisqr_diff)
        contour_chisqr_case = any(chisqr_cell_min .< chisqr_min .+ chisqr_contours .< chisqr_cell_max)
        return max_chisqr_case && diff_chisqr_case || contour_chisqr_case
    end

    function contour_cell_selector(i_cell::Int64, j_cell::Int64, grid::SimpleGrid)
        chisqr_min = grid.params[:chisqr_min]
        chisqr_cell = @view grid.value[:chisqr][i_cell:i_cell+1,j_cell:j_cell+1]
        chisqr_cell_min = minimum(chisqr_cell)
        chisqr_cell_max = maximum(chisqr_cell)
#        min_chisqr_case = chisqr_cell_min <= chisqr_min < chisqr_cell_max
        contour_chisqr_case = any(chisqr_cell_min .< chisqr_min .+ chisqr_contours .< chisqr_cell_max)
        return contour_chisqr_case
    end

    if tf.gsets.refinement_type == "nice"
        sell_selector = niceplot_cell_selector
    elseif tf.gsets.refinement_type == "contour"
        sell_selector = contour_cell_selector
    end

    function calculate_params!(grid::SimpleGrid)
        if haskey(grid.params, :chisqr_min)
            grid.params[:chisqr_min] = min(minimum(grid.value[:chisqr]), grid.params[:chisqr_min])
        else
            grid.params[:chisqr_min] = minimum(grid.value[:chisqr])
        end
        return nothing
    end

    if add_refinement == 0
        precalculate_Grid(tf.grid, get_ddstg_values_local, calculate_params!)
        for i in 1:tf.gsets.N_refinement
            tf.grid = refine_Grid(tf.grid, get_ddstg_values_local, sell_selector, calculate_params!)
        end
    else
        for i in 1:add_refinement
            tf.grid = refine_Grid(tf.grid, get_ddstg_values_local, sell_selector, calculate_params!)
        end
        tf.gsets = GridSetttings(tf.gsets.N_refinement + add_refinement, tf.gsets.CL, tf.gsets.contours, tf.gsets.refinement_type, tf.gsets.delta_chisqr_max, tf.gsets.delta_chisqr_diff, tf.gsets.gr_in_chisqr)
    end

    cut_ddstg_grid!(tf.grid)

    return tf
end