


struct EOSAgnosticTest{T}
    eos_list::Vector{String}
    tf_array::Vector{T} #tf - test framework
end

function EOSAgnosticTest(eos_list::Vector{String}, tf_sample::T) where {T}
    tf_array = [deepcopy(tf_sample) for i in 1:length(eos_list)]
    for i in 1:length(eos_list)
        tf_array[i].test = GeneralTest(tf_array[i].test, eos_list[i])
    end
    return(EOSAgnosticTest(eos_list,tf_array))
end

function Base.show(io::IO, eos_agn_test::EOSAgnosticTest)
    println(io, "Equation of state (EOS) agnostic test:")
    println(io,     "   List of EOS:  ", eos_agn_test.eos_list)
    print(io,       "   Sample Test:\n", eos_agn_test.tf_array[1])
	return nothing
end

function calculate!(eos_agn_test::EOSAgnosticTest; add_refinement=0)
#    N_refinement = eos_agn_test.tf_array[1].gsets.N_refinement
    work_dir = pwd()
    for i in 2:nworkers()+1
        rm("./worker$i", force=true, recursive=true)
        mkdir("./worker$i")
        cp(eos_agn_test.tf_array[1].tsets.par_file_init, "./worker$i/" * eos_agn_test.tf_array[1].tsets.par_file_init, force=true)
        cp(eos_agn_test.tf_array[1].tsets.tim_file, "./worker$i/" * eos_agn_test.tf_array[1].tsets.tim_file, force=true)
    end
    function single_task(tf::TempoFramework)
        cd(work_dir * "/worker$(myid())")
        calculate_t2!(tf, add_refinement=1)
    end
    for i in 0:(add_refinement != 0 ? add_refinement-1 : eos_agn_test.tf_array[1].gsets.N_refinement)
        eos_agn_test.tf_array .= pmap(single_task, eos_agn_test.tf_array)

        eos_agn_chisqr_grid = deepcopy(eos_agn_test.tf_array[1].grid.value[:chisqr])
        for j in 2:length(eos_agn_test.tf_array)
            eos_agn_chisqr_grid = min.(eos_agn_chisqr_grid, eos_agn_test.tf_array[j].grid.value[:chisqr])
        end
        for j in 1:length(eos_agn_test.tf_array)
            eos_agn_test.tf_array[j].grid.value[:eos_agn_chisqr] = eos_agn_chisqr_grid
            eos_agn_test.tf_array[j].grid.params[:eos_agn_chisqr_min] = minimum(eos_agn_chisqr_grid)
        end
        plot_test(eos_agn_test)
        savefig("eos_agn_test_$(maximum(eos_agn_test.tf_array[1].grid.ref_level)).png")
    end
end

function plot_test(eos_agn_test::EOSAgnosticTest)
    #    for i in 1:length(eos_agn_test.tf_array)
    #        cut_ddstg_grid!(eos_agn_test.tf_array[i].grid)
    #    end
    
        tf = eos_agn_test.tf_array[1]
    
        cm = ColorMap(ColorSchemes.okabe_ito.colors, 8)
        rc("mathtext",fontset="cm")
        rc("font", family="serif", size=12)
        fig, ax = subplots()eo
        pclm = ax.pcolormesh(tf.grid.y, tf.grid.x, tf.grid.value[:eos_agn_chisqr] .- tf.grid.params[:eos_agn_chisqr_min], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=tf.gsets.delta_chisqr_max), rasterized=true)
    
        for i in 1:length(eos_agn_test.tf_array)
            tf = eos_agn_test.tf_array[i]
            cs = ax.contour(tf.grid.y, tf.grid.x, round.(tf.grid.value[:chisqr] .- tf.grid.params[:chisqr_min], digits=1), levels=tf.gsets.contours, colors=[cm(i)])
            #clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
            plot([], [], label=tf.test.eosname, color=cm(i))
        end
    #    ax.set_ylabel(get_label("log10alpha0"), size=16)
    #    ax.set_xlabel(get_label("beta0"), size=16)
    
        legend(fontsize=11)
        cbar = colorbar(pclm)
    #    cbar.set_label(L"$\Delta\chi^{2}$", size=14)
        #title(L"0.1 \mu s")
        #ax.invert_xaxis()
        tight_layout()
    end
    