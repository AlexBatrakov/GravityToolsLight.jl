


struct EOSAgnosticTest{T}
    eos_list::Vector{String}
    tf_array::Vector{T} #tf - test framework
end

function EOSAgnosticTest(eos_list::Vector{String}, tf_sample::T) where {T}
    tf_array = [copy(tf_sample) for i in 1:length(eos_list)]
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

function calculate!(eos_agn_test::EOSAgnosticTest)
    work_dir = pwd()
    for i in 2:nworkers()+1
        mkdir("worker$i")
        mv(eos_agn_test.tf_array[1].tsets.par_file_init, "./"worker$i"/" * eos_agn_test.tf_array[1].tsets.par_file_init, force=true)
        mv(eos_agn_test.tf_array[1].tsets.tim_file, "./"worker$i"/" * eos_agn_test.tf_array[1].tsets.tim_file, force=true)
    end
end
