


struct EOSAgnosticTest{T}
    eos_list::Vector{String}
    tf::T #tf - test framework
end

calculate!(eos_agn_test::EOSAgnosticTest)
    tf_array = Array{typeof(eos_agn_test.tf)}(undef, length(eos_agn_test.eos_list))
    for i in 1:length(eos_agn_test.eos_list)
        tf_array[i].test.eosname = eos_agn_test.eos_list[i]
    end
end