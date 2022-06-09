struct MMDiagram
	m1::Vector{Float64}
	m2::NamedTuple{(:u, :m, :d), Tuple{StructArray{PKType},StructArray{PKType},StructArray{PKType}}}
    K_params_obs::KObsType
    PK_params_obs::PKObsType
	function MMDiagram(m1, K_params_obs, PK_params_obs)
		m2 = (u = StructArray{PKType}(undef, 0), m = StructArray{PKType}(undef, 0), d = StructArray{PKType}(undef, 0))
		return new(m1, m2, K_params_obs, PK_params_obs)
	end
end

function calculate!(mm::MMDiagram, pf::PhysicalFramework)
    function find_comp_mass!(F, x, m1, PK_name, PK_obs_value)
        pf.bnsys.psr.mass = m1
        pf.bnsys.comp.mass = abs(x[1])
        interpolate_bnsys!(pf)
        F[1] = (pf.bnsys.PK_params[PK_name] / PK_obs_value) - 1.0
    end

    for PK_name in PK_list
        for (i_m1, m1) in enumerate(mm.m1)
            println("pulsar mass = $m1, PK = $(PK_name)")
            PK_val, PK_err = mm.PK_params_obs[PK_name].val, mm.PK_params_obs[PK_name].err
            PK_obs_values = (PK_val + PK_err, PK_val, PK_val - PK_err)
            for (i_line, PK_obs_value) in enumerate(PK_obs_values)
                find_comp_mass_local!(F, x) = find_comp_mass!(F, x, m1, PK_name, PK_obs_value)
                solution = nlsolve(find_comp_mass_local!, [pf.bnsys.comp.mass])
                m2 = solution.zero[1]
                pf.bnsys.comp.mass = m2
                mm.m2[i_line][PK_name][i_m1] = m2
            end
        end
    end
    return mm, pf
end