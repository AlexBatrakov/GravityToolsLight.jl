KObsType = NamedTuple{K_list, NTuple{length(K_list), Measurement{Float64}}}
PKObsType = NamedTuple{PK_list, NTuple{length(PK_list), Measurement{Float64}}}
XObsType = NamedTuple{X_list, NTuple{length(X_list), Measurement{Float64}}}

mutable struct ClassicalTest
    name::String
    K_params_obs::KObsType
    PK_params_obs::PKObsType
    X_params_obs::XObsType
    N_refinement::Int64
    CL::Float64
    grid::SimpleGrid
    function ClassicalTest(name, K_params_obs, PK_params_obs)
        return new(name, K_params_obs, PK_params_obs, XObsType(zeros(length(X_list))), 6, 0.90)
    end
    function ClassicalTest(name, K_params_obs, PK_params_obs, X_params_obs)
        return new(name, K_params_obs, PK_params_obs, X_params_obs, 6, 0.90)
    end
end

function Base.show(io::IO, ct::ClassicalTest)
    println(io, "ClassicalTest")
    println(io, "Name of the pulsar: ", ct.name)
	println(io, "Observed Keplerian parameters: ", ct.K_params_obs)
	println(io, "Observed Post-Keplerian parameters: ", ct.PK_params_obs)
	println(io, "Observed Extra parameters: ", ct.X_params_obs)
	println(io, "Desired refinement level: ", ct.N_refinement)
    println(io, "Desired confidence level: ", ct.CL)
	return nothing
end

function find_initial_masses(ct::ClassicalTest, pf::DEFPhysicalFramework)
    PK_first = :k
    PK_second = :gamma
    for PK in keys(ct.PK_params_obs)
        if abs(ct.PK_params_obs[PK].val / ct.PK_params_obs[PK].err) > abs(ct.PK_params_obs[PK_first].val / ct.PK_params_obs[PK_first].err)
            PK_second = PK_first
            PK_first = PK
        elseif abs(ct.PK_params_obs[PK].val / ct.PK_params_obs[PK].err) > abs(ct.PK_params_obs[PK_second].val / ct.PK_params_obs[PK_second].err) && PK != PK_first
            PK_second = PK
        end
    end

    println(PK_first, " ", PK_second)

    function find_intersection!(F, x)
        pf.bnsys.psr.mass = abs(x[1])
        pf.bnsys.comp.mass = abs(x[2])
        interpolate_bnsys!(pf)
        F[1] = (pf.bnsys.PK_params[PK_first] / ct.PK_params_obs[PK_first].val) - 1.0 #/ ct.PK_params_sigma[PK_first]
        F[2] = (pf.bnsys.PK_params[PK_second] / ct.PK_params_obs[PK_second].val) - 1.0 #/ ct.PK_params_sigma[PK_first]
    end
    solution = nlsolve(find_intersection!, [pf.bnsys.psr.mass; pf.bnsys.comp.mass])
#    println(solution)
    m1, m2 = solution.zero
    pf.bnsys.psr.mass = m1
    pf.bnsys.comp.mass = m2
    return (m1 = m1, m2 = m2)
end

function get_chisqr(ct::ClassicalTest, pf::DEFPhysicalFramework)
    chisqr = 0.0
    for PK in keys(ct.PK_params_obs)
        if ct.PK_params_obs[PK].err != 0.0
            chisqr += ((pf.bnsys.PK_params[PK] - ct.PK_params_obs[PK].val) / ct.PK_params_obs[PK].err)^2
        end
    end
    for X in keys(ct.X_params_obs)
        if ct.X_params_obs[X].err != 0.0
            chisqr += ((pf.bnsys.X_params[X] - ct.X_params_obs[X].val) / ct.X_params_obs[X].err)^2
        end
    end
    return chisqr
end

function check_terms_in_chisqr(ct::ClassicalTest, pf::DEFPhysicalFramework)
    println("(theory - observation) / sigma")
    for PK in keys(ct.PK_params_obs)
        if ct.PK_params_obs[PK].err != 0.0
            println("$PK: ", ((pf.bnsys.PK_params[PK] - ct.PK_params_obs[PK].val) / ct.PK_params_obs[PK].err))
        end
    end
    for X in keys(ct.X_params_obs)
        if ct.X_params_obs[X].err != 0.0
            println("$X: ", ((pf.bnsys.X_params[X] - ct.X_params_obs[X].val) / ct.X_params_obs[X].err))
        end
    end
end

function find_best_masses(ct::ClassicalTest, pf::DEFPhysicalFramework)
    function get_chisqr_local(x)
        pf.bnsys.psr.mass = abs(x[1])
        pf.bnsys.comp.mass = abs(x[2])
        interpolate_bnsys!(pf)
        return get_chisqr(ct::ClassicalTest, pf::DEFPhysicalFramework)
    end
    find_initial_masses(ct, pf)
    sol = optimize(get_chisqr_local, [pf.bnsys.psr.mass, pf.bnsys.comp.mass])
    m1, m2 = Optim.minimizer(sol)
    pf.bnsys.psr.mass, pf.bnsys.comp.mass = m1, m2
    interpolate_bnsys!(pf)
    return (m1 = m1, m2 = m2)
end


function calculate!(ct::ClassicalTest, pf::DEFPhysicalFramework)

    function get_chisqr_local(log10alpha0::Float64, beta0::Float64)
        pf.theory.alpha0 = -exp10(log10alpha0)
        pf.theory.beta0  = beta0
        interpolate_mgrid!(pf)
        find_best_masses(ct, pf)
        chisqr = get_chisqr(ct, pf)
        println("alpha0 = $(pf.theory.alpha0), beta0 = $(pf.theory.beta0), m1 = $(pf.bnsys.psr.mass), m2 = $(pf.bnsys.comp.mass), chisqr = $chisqr")
        return ((:chisqr, :m1, :m2, :Pbdot),  (chisqr, pf.bnsys.psr.mass, pf.bnsys.comp.mass, pf.bnsys.PK_params.Pbdot) )
    end

    lvl = quantile(Chisq(2), ct.CL)

    function niceplot_cell_selector(i_cell::Int64, j_cell::Int64, grid::SimpleGrid)
        lvl = quantile(Chisq(2), 0.90)
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

    precalculate_Grid(ct.grid, get_chisqr_local, calculate_params!)
    for i in 1:ct.N_refinement
        ct.grid = refine_Grid(ct.grid, get_chisqr_local, niceplot_cell_selector, calculate_params!)
    end

    return ct, pf
end

ct_dataset = Dict{String, ClassicalTest}()

ct_dataset["Double Pulsar"] = ClassicalTest(
    "J0737−3039A",
    (Pb = 0.1022515592973 ± 0.0000000000010,
    T0= 55700.233017540 ± 0.000000013,
    e0 = 0.087777023 ± 0.000000061,
    omega0 = 204.753686 ± 0.000047,
    x0 = 1.415028603 ± 0.000000092),
    (k = (16.899323 ± 0.000013) / 360 * 0.1022515592973/365.25,
    gamma = 0.384045e-3 ± 0.000094e-3,
    Pbdot = -1.247752e-12 ± 0.000079e-12,
    r = 6.162e-6 ± 0.021e-6,
    s = 0.999936 ± 0.000010))

ct_dataset["J2222−0137"] = ClassicalTest(
    "J2222−0137",
    (Pb = 2.445759995471 ± 0.000000000006,
    T0= 58001.200912600 ± 0.000000002,
    e0 = 0.00032836 ± 0.00000002,
    omega0 = 120.45820 ± 0.00143,
    x0 = 10.8480213 ± 0.0000002),
    (k = (0.09607 ± 0.00048) / 360 * 2.445759995471/365.25,
    gamma = 0.0 ± 0.0,
#   Pbdot = -0.0143e-12 ± 0.0076e-12,
    Pbdot = -0.0080e-12 ± 0.0076e-12,
    r = (1.312 ± 0.009)*G*M_sun/c^3,
    s = 0.996623 ± 0.000115))