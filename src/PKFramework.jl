#K_list = (:Pb, :T0, :e0, :omega0, :x0)
#PK_list = (:k, :gamma, :Pbdot, :r, :s, :h3, :varsigma, :dtheta)
#X_list = (:m2, :q, :deltaN)

KObsType = NamedTuple{K_list, NTuple{length(K_list), Measurement{Float64}}}
PKObsType = NamedTuple{PK_list, NTuple{length(PK_list), Measurement{Float64}}}
XObsType = NamedTuple{X_list, NTuple{length(X_list), Measurement{Float64}}}

struct ObsParams
    K::KObsType
    PK::PKObsType
    X::XObsType
end

function Base.show(io::IO, params::ObsParams)
    println(io, "Observed parameters:")
    println(io, "Observed Keplerian parameters:\n   ", params.K)
	println(io, "Observed Post-Keplerian parameters:\n   ", params.PK)
	print(io,   "Observed Extra parameters:\n   ", params.X)
	return nothing
end


mutable struct PKFramework{T <: AbstractGravityTest}
    test::T
    obs_params::ObsParams
    gsets::GridSetttings
    grid::SimpleGrid
end

function Base.show(io::IO, pkf::PKFramework)
    println(io, "Post-Keplerian framework:")
    println(io, pkf.test)
    println(io, pkf.obs_params)
    print(io, pkf.gsets)
	return nothing
end


#PKFramework(;test::T, obs_params::ObsParams, gsets::GridSetttings) where {T <: AbstractGravityTest} = PKFramework{T}(test, obs_params, gsets, grid)

function PKFramework(test::GeneralTest, obs_params::ObsParams, gsets::GridSetttings)
    param1_grid = collect(LinRange(test.param1.min, test.param1.max, test.param1.N))
    param2_grid = collect(LinRange(test.param2.min, test.param2.max, test.param2.N))
    grid = SimpleGrid(Dict(), param1_grid, param2_grid)
    return PKFramework(test, obs_params, gsets, grid)
end

function find_initial_masses(obs_params::ObsParams, pf::DEFPhysicalFramework)
    PK_first = :k
    PK_second = :gamma
    for PK in keys(obs_params.PK)
        if abs(obs_params.PK[PK].val / obs_params.PK[PK].err) > abs(obs_params.PK[PK_first].val / obs_params.PK[PK_first].err)
            PK_second = PK_first
            PK_first = PK
        elseif abs(obs_params.PK[PK].val / obs_params.PK[PK].err) > abs(obs_params.PK[PK_second].val / obs_params.PK[PK_second].err) && PK != PK_first
            PK_second = PK
        end
    end

#    println(PK_first, " ", PK_second)

    function find_intersection!(F, x)
        pf.bnsys.psr.mass = abs(x[1])
        pf.bnsys.comp.mass = abs(x[2])
        interpolate_bnsys!(pf)
        F[1] = (pf.bnsys.PK_params[PK_first] / obs_params.PK[PK_first].val) - 1.0 #/ ct.PK_params_sigma[PK_first]
        F[2] = (pf.bnsys.PK_params[PK_second] / obs_params.PK[PK_second].val) - 1.0 #/ ct.PK_params_sigma[PK_first]
    end
    solution = nlsolve(find_intersection!, [pf.bnsys.psr.mass; pf.bnsys.comp.mass])
#    println(solution)
    m1, m2 = abs.(solution.zero)
    pf.bnsys.psr.mass = m1
    pf.bnsys.comp.mass = m2
    return (m1 = m1, m2 = m2)
end

function get_chisqr(obs_params::ObsParams, pf::DEFPhysicalFramework)
    chisqr = 0.0
    for PK in keys(obs_params.PK)
        if obs_params.PK[PK].err != 0.0
            chisqr += ((pf.bnsys.PK_params[PK] - obs_params.PK[PK].val) / obs_params.PK[PK].err)^2
        end
    end
    for K in keys(obs_params.K)
        if obs_params.K[K].err != 0.0
            chisqr += ((pf.bnsys.K_params[K] - obs_params.K[K].val) / obs_params.K[K].err)^2
        end
    end
    for X in keys(obs_params.X)
        if obs_params.X[X].err != 0.0
            chisqr += ((pf.bnsys.X_params[X] - obs_params.X[X].val) / obs_params.X[X].err)^2
        end
    end
    return chisqr
end

function check_terms_in_chisqr(obs_params::ObsParams, pf::DEFPhysicalFramework)
    println("(theory - observation) / sigma")
    for PK in keys(obs_params.PK)
        if obs_params.PK[PK].err != 0.0
            println("$PK: ", ((pf.bnsys.PK_params[PK] - obs_params.PK[PK].val) / obs_params.PK[PK].err))
        end
    end
    for K in keys(obs_params.K)
        if obs_params.K[K].err != 0.0
            println("$K: ", ((pf.bnsys.K_params[K] - obs_params.K[K].val) / obs_params.K[K].err))
        end
    end
    for X in keys(obs_params.X)
        if obs_params.X[X].err != 0.0
            println("$X: ", ((pf.bnsys.X_params[X] - obs_params.X[X].val) / obs_params.X[X].err))
        end
    end
end

function optimize_PK_method(obs_params::ObsParams, pf::DEFPhysicalFramework)
    function get_chisqr_local(x)
        pf.bnsys.psr.mass = abs(x[1])
        pf.bnsys.comp.mass = abs(x[2])
        Pb = x[3]
        e0 = x[4]
        x0 = x[5]
        pf.bnsys.K_params = (Pb = Pb, T0 = pf.bnsys.K_params.T0, e0 = e0, omega0 = pf.bnsys.K_params.omega0, x0 = x0)
        interpolate_bnsys!(pf)
        return sqrt(get_chisqr(obs_params, pf))
    end
    #find_initial_masses(obs_params, pf)
    sol = optimize(get_chisqr_local, [pf.bnsys.psr.mass, pf.bnsys.comp.mass, pf.bnsys.K_params.Pb, pf.bnsys.K_params.e0, pf.bnsys.K_params.x0])
    m1, m2, Pb, e0, x0 = Optim.minimizer(sol)
    m1 = abs(m1); m2 = abs(m2)
    pf.bnsys.psr.mass, pf.bnsys.comp.mass = m1, m2
    pf.bnsys.K_params = (Pb = Pb, T0 = pf.bnsys.K_params.T0, e0 = e0, omega0 = pf.bnsys.K_params.omega0, x0 = x0)
    interpolate_bnsys!(pf)
    return (m1 = m1, m2 = m2, Pb = Pb, e0 = e0,  x0 = x0)
end

function find_best_masses(obs_params::ObsParams, pf::DEFPhysicalFramework)
    function get_chisqr_local(x)
        pf.bnsys.psr.mass = abs(x[1])
        pf.bnsys.comp.mass = abs(x[2])
        interpolate_bnsys!(pf)
        return sqrt(get_chisqr(obs_params, pf))
    end
    #find_initial_masses(obs_params, pf)
    sol = optimize(get_chisqr_local, [pf.bnsys.psr.mass, pf.bnsys.comp.mass])
    m1, m2 = abs.(Optim.minimizer(sol))
    pf.bnsys.psr.mass, pf.bnsys.comp.mass = m1, m2
    interpolate_bnsys!(pf)
    return (m1 = m1, m2 = m2)
end

#function get_chisqr(fixed_params::Dict, pf::DEFPhysicalFramework)
#    chisqr = 0.0
#    for name in keys(fixed_params)
#        if name ∈ keys(pf.bnsys.PK_params)
#            chisqr += ((pf.bnsys.PK_params[name] - fixed_params[name].val) / fixed_params[name].err)^2
#        elseif name ∈ keys(pf.bnsys.X_params)
#            chisqr += ((pf.bnsys.X_params[name] - fixed_params[name].val) / fixed_params[name].err)^2
#        elseif name == :m1
#            chisqr += ((pf.bnsys.psr.mass -  - fixed_params[name].val) / fixed_params[name].err)^2
#        elseif name == :m2
#            chisqr += ((pf.bnsys.comp.mass -  - fixed_params[name].val) / fixed_params[name].err)^2
#        end
#    end
#    return chisqr
#end
#
#function find_best_fit(fixed_params::Dict, pf::DEFPhysicalFramework)
#    function get_chisqr_local(x)
#        pf.bnsys.psr.mass = abs(x[1])
#        pf.bnsys.comp.mass = abs(x[2])
#        interpolate_bnsys!(pf)
#        chisqr = get_chisqr(fixed_params, pf)
#        return sqrt(chisqr)
#    end
#    sol = optimize(get_chisqr_local, [pf.bnsys.psr.mass, pf.bnsys.comp.mass])
#    println(sol)
#    m1, m2 = abs.(Optim.minimizer(sol))
#    pf.bnsys.psr.mass, pf.bnsys.comp.mass = m1, m2
#    interpolate_bnsys!(pf)
#    return (m1 = m1, m2 = m2)
#end

function find_masses(params, pf::DEFPhysicalFramework)
    name1, name2 = keys(params)
    function get_residuals!(F, x)
        pf.bnsys.psr.mass = abs(x[1])
        pf.bnsys.comp.mass = abs(x[2])
        interpolate_bnsys!(pf)
        F[1] = (pf.bnsys.PK_params[Symbol(name1)] - params[name1]) / (params[name1] == 0 ? 1 : params[name1])
        F[2] = (pf.bnsys.PK_params[Symbol(name2)] - params[name2]) / (params[name2] == 0 ? 1 : params[name2])
    end
    sol = nlsolve(get_residuals!, [pf.bnsys.psr.mass, pf.bnsys.comp.mass])
    m1, m2 = abs.(sol.zero)
    pf.bnsys.psr.mass, pf.bnsys.comp.mass = m1, m2
    interpolate_bnsys!(pf)
    if !converged(sol)
        m1, m2 = NaN, NaN
    end
    return (m1 = m1, m2 = m2)
end

function calculate!(pkf::PKFramework, pf::DEFPhysicalFramework)

    pf.theory.alpha0 = 0.0
    pf.theory.beta0  = 0.0
    m1_GR, m2_GR = find_best_masses(pkf.obs_params, pf)

    if typeof(pkf.test.alpha0) == Float64
        pf.theory.alpha0 = pkf.test.alpha0
    end
    if typeof(pkf.test.beta0) == Float64
        pf.theory.beta0 = pkf.test.beta0
    end
    interpolate_mgrid!(pf)

    key_change_gravity = pkf.test.param1.name == "log10alpha0" || pkf.test.param1.name == "alpha0" || pkf.test.param1.name == "beta0"

    pkf.grid.params[:chisqr_min] = Inf

    params = Dict(pkf.test.param1.name => pkf.test.param1.min, pkf.test.param2.name => pkf.test.param2.min)

    function get_chisqr_local(param1::Float64, param2::Float64)
        name1, name2 = pkf.test.param1.name, pkf.test.param2.name
        key_fit_masses = false
        if key_change_gravity == true
            key_fit_masses = true
            if name1 == "log10alpha0"
                pf.theory.alpha0 = -exp10(param1)
            elseif  name1 == "alpha0"
                pf.theory.alpha0 = param1
            elseif name1 == "beta0"
                pf.theory.beta0  = param1
            end
            if name2 == "log10alpha0"
                pf.theory.alpha0 = -exp10(param2)
            elseif  name2 == "alpha0"
                pf.theory.alpha0 = param2
            elseif name2 == "beta0"
                pf.theory.beta0  = param2
            end
            interpolate_mgrid!(pf)
        end
        pf.bnsys.psr.mass = m1_GR
        pf.bnsys.comp.mass = m2_GR
        if key_fit_masses == true
            interpolate_mgrid!(pf)
            find_best_masses(pkf.obs_params, pf)
        else
            find_masses((name1 = param1, name2 = param2), pf)
        end
        chisqr = get_chisqr(pkf.obs_params, pf)
        @printf "run %s = %10.6f, %s = %10.6f, χ2 = %10.3f\n" pkf.test.param1.name param1 pkf.test.param2.name param2 chisqr
        return ((:chisqr, :m1, :m2, :Pbdot, :alphaA, :betaA, :kA),  (chisqr, pf.bnsys.psr.mass, pf.bnsys.comp.mass, pf.bnsys.PK_params.Pbdot, pf.bnsys.psr.alphaA, pf.bnsys.psr.betaA, pf.bnsys.psr.kA) )
    end

    lvl = quantile(Chisq(2), pkf.gsets.CL)

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

    precalculate_Grid(pkf.grid, get_chisqr_local, calculate_params!)
    for i in 1:pkf.gsets.N_refinement
        pkf.grid = refine_Grid(pkf.grid, get_chisqr_local, niceplot_cell_selector, calculate_params!)
    end

    return pkf, pf
end


obs_params_dataset = Dict{String, ObsParams}()

obs_params_dataset["J1141-6545_DD"] = ObsParams(
    (Pb = 0.19765095844298 ± 0.00000000003768,
    T0= 53999.9958600347 ± 0.0000002756,
    e0 = 0.1718768077 ± 0.0000030767,
    omega0 = 80.691561013198 ± 0.002501062255,
    x0 = 1.858915489 ± 0.000016012),
    (k = (5.3105030 ± 0.0001169) / 360 * 0.19765095844298/365.25,
    gamma = 0.000783621 ± 0.000080708,
#    Pbdot = -0.3855839e-12 ± 0.0168915e-12,
    Pbdot = -0.3863656e-12 ± 0.0168915e-12,
    r = (0.866092 ± 0.752243)*G*M_sun/c^3,
    s = 0.970686 ± 0.057246,
    h3 = 0.0 ± 0.0,
    varsigma = 0.0 ± 0.0,
    dtheta = 0.0 ± 0.0),
    (m2 = 0.0 ± 0.0,
    q = 0.0 ± 0.0,
    deltaN = 0.0)
    )

    obs_params_dataset["J1952+2630_DDFWHE_old"] = ObsParams(
    (Pb = 0.39188059665094 ± 0.00000029819894,
    T0= 55407.6179751086 ± 0.0022389933,
    e0 = 0.0000413989 ± 0.0000013364,
    omega0 = 289.735318984806 ± 2.056801088852,
    x0 = 2.798179920 ± 0.000007079),
    (k = (1.6814251 ± 0.2553355) / 360 * 0.39188059665094/365.25,
    gamma = 0.000001764 ± 0.000011182,
#    Pbdot = -0.0997964e-12 ± 0.0372482e-12,
    Pbdot = -0.0997964e-12 + 0.0018e-12 ± 0.0372482e-12,
    r = (0.0 ± 0.0)*G*M_sun/c^3,
    s = 0.0 ± 0.0,
    h3 = 0.000001306850 ± 0.000000589634,
    varsigma = 0.0 ± 0.0,
    dtheta = 0.0 ± 0.0),
    (m2 = 0.0 ± 0.0,
    q = 0.0 ± 0.0,
    deltaN = 0.0)
    )

obs_params_dataset["J1952+2630_DDFWHE"] = ObsParams(
    (Pb = 0.39188062491079 ± 0.00000030762645,
    T0= 55407.6182002729 ± 0.0022693526,
    e0 = 0.0000409978 ± 0.0000013932,
    omega0 = 289.942164350072 ± 2.084690036136,
    x0 = 2.798183184 ± 0.000007391),
    (k = (1.7056226 ± 0.2634040) / 360 * 0.39188059665094/365.25,
    gamma = 0.000001878 ± 0.000011185,
#    Pbdot = -0.0992974e-12 ± 0.0463139e-12,
    Pbdot = -0.0992974e-12 + 0.0018e-12 ± 0.0463139e-12,
    r = (0.0 ± 0.0)*G*M_sun/c^3,
    s = 0.0 ± 0.0,
    h3 = 0.000001308320 ± 0.000000592080,
    varsigma = 0.0 ± 0.0,
    dtheta = 0.0 ± 0.0),
    (m2 = 0.0 ± 0.0,
    q = 0.0 ± 0.0,
    deltaN = 0.0)
    )

