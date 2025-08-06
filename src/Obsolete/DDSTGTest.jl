struct DDSTGTestSettings
    par_file_init::String
    tim_file::String
    add_flag::String
    fit_XPBDOT::Bool
end

DDSTGTestSettings(;par_file_init, tim_file, add_flag, fit_XPBDOT) = DDSTGTestSettings(par_file_init, tim_file, add_flag, fit_XPBDOT)

mutable struct DDSTGTest
    psrname::String
    eosname::String
    N_refinement::Int64
    CL::Float64
    sets::DDSTGTestSettings
    grid::SimpleGrid
    function DDSTGTest(psrname, eosname, N_refinement, CL, sets)
        return new(psrname, eosname, N_refinement, CL, sets)
    end
end

DDSTGTest(;psrname, eosname, N_refinement, CL, sets) = DDSTGTest(psrname, eosname, N_refinement, CL, sets)


function Base.show(io::IO, ddstgt::DDSTGTest)
    println(io, "ClassicalTest")
    println(io, "Name of the pulsar: ", ddstgt.psrname)
    println(io, "Equation of state: ", ddstgt.eosname)
	println(io, "Desired refinement level: ", ddstgt.N_refinement)
    println(io, "Desired confidence level: ", ddstgt.CL)
	return nothing
end

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
    return chisqr
end

function read_params(params, par_file)
    open(par_file,"r") do file_out
        for line in eachline(par_file)
            for param_name in keys(params)
                if startswith(line, String(param_name)*" ")
                    try
                        params[param_name] = parse(Float64,line[11:26])
#                    println(line)
                    catch err
                        params[param_name] = NaN
                    end
                end
            end
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
    write_new_params((NITS=30, XPBDOT=fit_XPBDOT), par_file_work)

    run_tempo(par_file_work, tim_file; silent=true, add_flag="-c "*add_flag)

    cp(par_file_out, par_file_work; force=true)
    write_new_params((NITS=5, XPBDOT=false), par_file_work)
    chisqr = read_chisqr_tempo_lis()
#    println(chisqr)
    return mean(chisqr[10:30])
end

function calculate!(ddstgt::DDSTGTest)

    par_file_init = ddstgt.sets.par_file_init
    par_file_work = "$(par_file_init[1:end-4])_work.par"
    par_file_out = "$(ddstgt.psrname).par"
    tim_file = ddstgt.sets.tim_file
    add_flag = ddstgt.sets.add_flag
    eosname = ddstgt.eosname
    fit_XPBDOT = ddstgt.sets.fit_XPBDOT
    lvl = quantile(Chisq(2), ddstgt.CL)

    println("Obtaining a working parfile for fit_XPBDOT=$fit_XPBDOT")
    ddstgt.grid.params[:chisqr_min] = get_par_file_work(par_file_init, par_file_out, tim_file; add_flag=add_flag, fit_XPBDOT=fit_XPBDOT)
    
    function get_ddstg_values_local(log10alpha0, beta0; silent=true)
        write_new_params((ALPHA0=-exp10(log10alpha0), BETA0=beta0, EOS=eosname, NITS=5, GAIN=1.0), par_file_work)
        chisqr = run_tempo(par_file_work, tim_file; silent=true, add_flag=add_flag)
        if isinf(chisqr)
            println("The first attempt failed. Reajustment of parameters.")
            modify_par_file(Dict("NITS"=>10, "GAIN"=>0.1), par_file_work)
            run_tempo(par_file_work, tim_file; silent=true, add_flag="-c "*add_flag)
            modify_par_file(Dict("NITS"=>5, " GAIN"=>1.0), par_file_out)
            chisqr = run_tempo(par_file_out, tim_file; silent=true, add_flag=add_flag)
        end
        @printf "run α0 = %10.6f, β0 = %10.6f, χ2 = %10.3f\n" -exp10(log10alpha0) beta0 chisqr
        temp_dict = read_params(Dict(:A1=>0.0, :E=>0.0, :T0=>0.0, :PB=>0.0, :OM=>0.0, :OMDOT=>0.0, :GAMMA=>0.0, :PBDOT=>0.0, :SINI=>0.0, :DTHETA=>0.0, :XDOT=>0.0, :DR=>0.0,:MA=>0.0, :MB =>0.0, :ALPHA0=>0.0, :BETA0=>0.0, :ALPHAA=>0.0, :BETAA=>0.0, :kA=>0.0), par_file_out)
        ddstg_names = tuple(:chisqr, keys(temp_dict)...)
        ddstg_values = tuple(chisqr, values(temp_dict)...)
        return (ddstg_names, ddstg_values)
    end

    function niceplot_cell_selector(i_cell::Int64, j_cell::Int64, grid::SimpleGrid)
        δχ2_max = lvl_5σ = quantile(Chisq(2), 0.999999426696856)
        δχ2_min = 0.5
        χ2_min = grid.params[:chisqr_min]
        χ2_cell = @view grid.value[:chisqr][i_cell:i_cell+1,j_cell:j_cell+1]
        χ2_cell_min = minimum(χ2_cell)
        χ2_cell_max = maximum(χ2_cell)
        large_χ2_case = (χ2_cell_min < χ2_min + δχ2_max)
        diff_χ2_case = (χ2_cell_max - χ2_cell_min > δχ2_min)
        contour_χ2_case = χ2_cell_min < χ2_min + lvl < χ2_cell_max
        return large_χ2_case && diff_χ2_case || contour_χ2_case
    end

    function calculate_params!(grid::SimpleGrid)
        if haskey(grid.params, :chisqr_min)
            grid.params[:chisqr_min] = min(minimum(grid.value[:chisqr]), grid.params[:chisqr_min])
        else
            grid.params[:chisqr_min] = minimum(grid.value[:chisqr])
        end
        return nothing
    end

    precalculate_Grid(ddstgt.grid, get_ddstg_values_local, calculate_params!)
    for i in 1:ddstgt.N_refinement
        ddstgt.grid = refine_Grid(ddstgt.grid, get_ddstg_values_local, niceplot_cell_selector, calculate_params!)
    end

    return ddstgt
end