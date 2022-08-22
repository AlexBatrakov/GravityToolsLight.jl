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
    masses_init::NamedTuple{(:m1, :m2), Tuple{Float64, Float64}}
end

ObsParams(;Pb = 0.0 ± 0.0, T0 = 0.0 ± 0.0, e0 = 0.0 ± 0.0, eps1 = 0.0 ± 0.0, eps2 = 0.0 ± 0.0, omega0 = 0.0 ± 0.0, x0 = 0.0 ± 0.0, k = 0.0 ± 0.0, omegadot = 0.0 ± 0.0, gamma = 0.0 ± 0.0, Pbdot = 0.0 ± 0.0, r = 0.0 ± 0.0, m2_shapiro = 0.0 ± 0.0, s = 0.0 ± 0.0, h3 = 0.0 ± 0.0, varsigma = 0.0 ± 0.0, dtheta = 0.0 ± 0.0, m2 = 0.0 ± 0.0, q = 0.0 ± 0.0, deltaN = 0.0 ± 0.0, m1_init = 1.0, m2_init = 1.0) = ObsParams(KObsType((Pb, T0, e0 == 0 ? sqrt(eps1^2 + eps2^2) : e0, omega0 == 0 ? 180/pi*atan(eps1/eps2) : omega0, x0)), PKObsType((k == 0 ? omegadot / 360 * Pb/365.25 : k, gamma, Pbdot, r == 0 ? m2_shapiro*G*M_sun/c^3 : r, s, h3, varsigma, dtheta)), XObsType((m2, q, deltaN)), (m1 = m1_init, m2 = m2_init))

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
    chisqr = isnan(chisqr) ? Inf : chisqr
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

function get_terms_in_chisqr(pkf::PKFramework, obs_params::ObsParams, pf::DEFPhysicalFramework)

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

function find_masses(params, obs_params, pf::DEFPhysicalFramework)
    name1, name2 = keys(params)
    function get_value(name)
        if name == "m1"
            value = pf.bnsys.psr.mass
        elseif name == "m2"
            value = pf.bnsys.comp.mass
        elseif haskey(pf.bnsys.PK_params, Symbol(name))
            value = pf.bnsys.PK_params[Symbol(name)]
        elseif haskey(pf.bnsys.X_params, Symbol(name))
            value = pf.bnsys.X_params[Symbol(name)]
        elseif name == "cosi"
            s = pf.bnsys.PK_params[:s]
            value = sqrt(1 - s^2 >= 0 ? 1 - s^2 : NaN)
        end
    end

    function get_residuals!(F, x)
        pf.bnsys.psr.mass = abs(x[1])
        pf.bnsys.comp.mass = abs(x[2])
        interpolate_bnsys!(pf)
        F[1] = (get_value(name1) - params[name1]) / (params[name1] == 0 ? 1 : params[name1])
        F[2] = (get_value(name2) - params[name2]) / (params[name2] == 0 ? 1 : params[name2])
    end
    sol = nlsolve(get_residuals!, [obs_params.masses_init.m1, obs_params.masses_init.m2])
    m1, m2 = abs.(sol.zero)
    if !converged(sol)
        m1, m2 = NaN, NaN
    end
    pf.bnsys.psr.mass, pf.bnsys.comp.mass = m1, m2
    interpolate_bnsys!(pf)
    return (m1 = m1, m2 = m2)
end

function get_value(name, pf)
    if name == "m1"
        value = pf.bnsys.psr.mass
    elseif name == "m2" || name == "M2"
        value = pf.bnsys.comp.mass
    elseif haskey(pf.bnsys.PK_params, Symbol(name))
        value = pf.bnsys.PK_params[Symbol(name)]
    elseif haskey(pf.bnsys.X_params, Symbol(name))
        value = pf.bnsys.X_params[Symbol(name)]
    elseif name == "cosi" || name == "COSI"
        s = pf.bnsys.PK_params[:s]
        value = sqrt(1 - s^2 >= 0 ? 1 - s^2 : 1.0)
    elseif name == "PBDOT"
        s = pf.bnsys.PK_params[:Pbdot]
    elseif name == "GAMMA"
        s = pf.bnsys.PK_params[:gamma]
    end
end

function adjust_m1(params, pf::DEFPhysicalFramework)
    name1, name2 = keys(params)
    pf.bnsys.comp.mass = params["m2"]
    name = name1 == "m2" ? name2 : name1
    function get_residuals!(F, x)
        pf.bnsys.psr.mass = abs(x[1])
        interpolate_bnsys!(pf)
        F[1] = (get_value(name, pf) - params[name]) / (params[name] == 0 ? 1 : params[name])
    end
    sol = nlsolve(get_residuals!, [pf.bnsys.psr.mass])
    m1 = abs.(sol.zero)[1]
    if !converged(sol)
        m1 = NaN
    end
    pf.bnsys.psr.mass = m1
    interpolate_bnsys!(pf)
    return pf
end

function adjust_m1m2(params, pf::DEFPhysicalFramework)
    name1, name2 = keys(params)
    function get_residuals!(F, x)
        pf.bnsys.psr.mass = abs(x[1])
        pf.bnsys.comp.mass = abs(x[2])
        interpolate_bnsys!(pf)
        F[1] = (get_value(name1, pf) - params[name1]) / (params[name1] == 0 ? 1 : params[name1])
        F[2] = (get_value(name2, pf) - params[name2]) / (params[name2] == 0 ? 1 : params[name2])
    end
    sol = nlsolve(get_residuals!, [pf.bnsys.psr.mass, pf.bnsys.comp.mass])
    m1, m2 = abs.(sol.zero)
    if !converged(sol)
        m1, m2 = NaN, NaN
    end
    pf.bnsys.psr.mass, pf.bnsys.comp.mass = m1, m2
    interpolate_bnsys!(pf)
    return pf
end

function adjust_masses(params, pf::DEFPhysicalFramework)
    name1, name2 = keys(params)

    if haskey(params, "m1") && haskey(params, "m2")
        pf.bnsys.psr.mass = params["m1"]
        pf.bnsys.comp.mass = params["m2"]
        interpolate_bnsys!(pf)
    elseif haskey(params, "m1")
        pf.bnsys.psr.mass = params["m1"]
        name = name1 == "m1" ? name2 : name1
        function get_residuals!(F, x)
            pf.bnsys.comp.mass = abs(x[1])
            interpolate_bnsys!(pf)
            F[1] = (get_value(name, pf) - params[name]) / (params[name] == 0 ? 1 : params[name])
        end
        sol = nlsolve(get_residuals!, [pf.bnsys.comp.mass])
        m2 = abs.(sol.zero)[1]
        if !converged(sol)
            m2 = NaN
        end
        pf.bnsys.comp.mass = m2
        interpolate_bnsys!(pf)
    elseif haskey(params, "m2")
        adjust_m1(params, pf)
    else
        adjust_m1m2(params, pf)
    end

    return pf
end

function calculate!(pkf::PKFramework, pf::DEFPhysicalFramework; add_refinement=0)

    pf.theory.alpha0 = 0.0
    pf.theory.beta0  = 0.0
    m1_GR, m2_GR = find_best_masses(pkf.obs_params, pf)
    pkf.grid.params[:chisqr_gr] = get_chisqr(pkf.obs_params, pf)



    if typeof(pkf.test.alpha0) == Float64
        pf.theory.alpha0 = pkf.test.alpha0
    end
    if typeof(pkf.test.beta0) == Float64
        pf.theory.beta0 = pkf.test.beta0
    end
    interpolate_mgrid!(pf)

    key_change_gravity = pkf.test.param1.name == "log10alpha0" || pkf.test.param1.name == "alpha0" || pkf.test.param1.name == "beta0"

    if !haskey(pkf.grid.params, :chisqr_min)
        pkf.grid.params[:chisqr_min] = pkf.gsets.gr_in_chisqr ? pkf.grid.params[:chisqr_gr] : Inf
    end


    params = Dict(pkf.test.param1.name => pkf.test.param1.min, pkf.test.param2.name => pkf.test.param2.min)

    function get_chisqr_local(param1::Float64, param2::Float64)
        name1, name2 = pkf.test.param1.name, pkf.test.param2.name
        params[name1] = param1
        params[name2] = param2
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
            adjust_masses(params, pf)
        end
        chisqr = get_chisqr(pkf.obs_params, pf)
        pkf.grid.params[:chisqr_min] = pkf.grid.params[:chisqr_min] < chisqr ? pkf.grid.params[:chisqr_min] : chisqr
        @printf "run %s = %10.6f, %s = %10.6f, χ2 = %10.3f\n" pkf.test.param1.name param1 pkf.test.param2.name param2 chisqr
        return ((:chisqr, :m1, :m2, :k, :gamma, :Pbdot, :r, :s, :h3, :varsigma, :dtheta,  :alphaA, :betaA, :kA),  (chisqr, pf.bnsys.psr.mass, pf.bnsys.comp.mass, pf.bnsys.PK_params.k, pf.bnsys.PK_params.gamma, pf.bnsys.PK_params.Pbdot, pf.bnsys.PK_params.r, pf.bnsys.PK_params.s, pf.bnsys.PK_params.h3, pf.bnsys.PK_params.varsigma, pf.bnsys.PK_params.dtheta, pf.bnsys.psr.alphaA, pf.bnsys.psr.betaA, pf.bnsys.psr.kA) )
    end


    chisqr_contours = pkf.gsets.contours
    delta_chisqr_max = pkf.gsets.delta_chisqr_max
    delta_chisqr_diff = pkf.gsets.delta_chisqr_diff

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

    function contour_min_cell_selector(i_cell::Int64, j_cell::Int64, grid::SimpleGrid)
        chisqr_min = grid.params[:chisqr_min]
        chisqr_cell = @view grid.value[:chisqr][i_cell:i_cell+1,j_cell:j_cell+1]
        chisqr_cell_min = minimum(chisqr_cell)
        chisqr_cell_max = maximum(chisqr_cell)
        min_chisqr_case = chisqr_cell_min <= chisqr_min < chisqr_cell_max
        contour_chisqr_case = any(chisqr_cell_min .< chisqr_min .+ chisqr_contours .< chisqr_cell_max)
        return contour_chisqr_case || min_chisqr_case
    end

    function massmass_cell_selector(i_cell::Int64, j_cell::Int64, grid::SimpleGrid)
        PK_contour_case = false
        for PK in keys(pkf.obs_params.PK)
            if pkf.obs_params.PK[PK].err != 0.0
                PK_cell = @view grid.value[PK][i_cell:i_cell+1,j_cell:j_cell+1]
                PK_values = [pkf.obs_params.PK[PK].val - pkf.obs_params.PK[PK].err, pkf.obs_params.PK[PK].val, pkf.obs_params.PK[PK].val + pkf.obs_params.PK[PK].err]
                PK_cell_min = minimum(x -> isnan(x) ? +Inf : x, PK_cell)
                PK_cell_max = maximum(x -> isnan(x) ? -Inf : x, PK_cell)
                PK_contour_case = PK_contour_case || any(x -> PK_cell_min < x < PK_cell_max, PK_values)
                if any(isnan, PK_cell)
                    PK_contour_case = PK_contour_case || any(x -> PK_cell_min < x, PK_values) || any(x -> x < PK_cell_max, PK_values)
                end
                #println("$PK $(PK_contour_case) $(PK_cell_min) $(pkf.obs_params.PK[PK].val) $(PK_cell_max)")
            end
        end
        return PK_contour_case || contour_cell_selector(i_cell, j_cell, grid)
    end

    if pkf.gsets.refinement_type == "nice"
        sell_selector = niceplot_cell_selector
    elseif pkf.gsets.refinement_type == "contour"
        sell_selector = contour_cell_selector
    elseif pkf.gsets.refinement_type == "contour+min"
        sell_selector = contour_min_cell_selector
    elseif pkf.gsets.refinement_type == "massmass"
        sell_selector = massmass_cell_selector
    end

    function calculate_params!(grid::SimpleGrid)
        if haskey(grid.params, :chisqr_min)
            grid.params[:chisqr_min] = min(minimum(grid.value[:chisqr]), grid.params[:chisqr_min])
        else
            grid.params[:chisqr_min] = minimum(grid.value[:chisqr])
        end
        for PK in keys(pkf.obs_params.PK)
            if pkf.obs_params.PK[PK].err != 0.0
                grid.params[PK] = minimum(abs.(grid.value[PK] .- pkf.obs_params.PK[PK].val))
            end
        end
        return nothing
    end

    if add_refinement == 0
        precalculate_Grid(pkf.grid, get_chisqr_local, calculate_params!)
        for i in 1:pkf.gsets.N_refinement
            pkf.grid = refine_Grid(pkf.grid, get_chisqr_local, sell_selector, calculate_params!)
        end
    else
        for i in 1:add_refinement
            pkf.grid = refine_Grid(pkf.grid, get_chisqr_local, sell_selector, calculate_params!)
        end
        pkf.gsets = GridSetttings(pkf.gsets.N_refinement + add_refinement, pkf.gsets.CL, pkf.gsets.contours, pkf.gsets.refinement_type, pkf.gsets.delta_chisqr_max, pkf.gsets.delta_chisqr_diff, pkf.gsets.gr_in_chisqr)
    end

    return pkf, pf
end

function parse_par_file(par_file::String)
    pf = open(par_file, "r")
    pf_lines = readlines(pf)
    for line in pf_lines
        if startswith(line, "PB")
        end
    end
end


obs_params_dataset = Dict{String, ObsParams}()

obs_params_dataset["J1141-6545_DD"] = ObsParams(
    Pb = 0.19765095844297 ± 0.00000000003768,
    T0 = 53999.9958600347 ± 0.0000002756,
    e0 = 0.1718768077 ± 0.0000030767,
    omega0 = 80.691561090562 ± 0.002501062308,
    x0 = 1.858915488 ± 0.000016012,
    omegadot = 5.3105030 ± 0.0001169,
    gamma = 0.000783623 ± 0.000080708,
#    Pbdot = -0.3842285e-12 - 0.0013547e-12 ± 0.0168915e-12,
    Pbdot = -0.3842285e-12 ± 0.0168915e-12,
    m2_shapiro = 0.866098 ± 0.752274,
    s = 0.970685 ± 0.057250,
    m1_init = 1.273322,
    m2_init = 1.016393
    )

obs_params_dataset["J1141-6545_DDFWHE"] = ObsParams(
    Pb = 0.19765095844297 ± 0.00000000003768,
    T0 = 53999.9958600347 ± 0.0000002756,
    e0 = 0.1718768077 ± 0.0000030767,
    omega0 = 80.691561112205 ± 0.002501062293,
    x0 = 1.858915488 ± 0.000016012,
    omegadot = 5.3105030 ± 0.0001169,
    gamma = 0.000783624 ± 0.000080708,
#    Pbdot = -0.3842283e-12 - 0.0013547e-12 ± 0.0168915e-12,
    Pbdot = -0.3842283e-12 ± 0.0168915e-12,
    varsigma = 0.782585 ± 0.192029,
    h3 = 0.000002044640 ± 0.000000587384,
    m1_init = 1.273322,
    m2_init = 1.016393
    )

obs_params_dataset["J1141-6545_GAMMA_PBDOT"] = ObsParams(
    Pb = 0.19765095844297 ± 0.00000000003768,
    T0 = 53999.9958600347 ± 0.0000002756,
    e0 = 0.1718768077 ± 0.0000030767,
    omega0 = 80.691561112205 ± 0.002501062293,
    x0 = 1.858915488 ± 0.000016012,
    gamma = 0.000783624 ± 0.000080708,
#    Pbdot = -0.3842283e-12 - 0.0013547e-12 ± 0.0168915e-12,
    Pbdot = -0.3842283e-12 ± 0.0168915e-12,
    m1_init = 1.273322,
    m2_init = 1.016393
    )


obs_params_dataset["J1952+2630_DD"] = ObsParams(
    Pb = 0.39188062188725 ± 0.00000030768055,
    T0 = 55407.6182451646 ± 0.0022682574,
    e0 = 0.0000409993 ± 0.0000013868,
    omega0 = 289.983405516698 ± 2.083683808190,
    x0 = 2.798183224 ± 0.000007366,
    omegadot = 1.7030337 ± 0.2634503,
#   gamma = GAMMA          0.000000225 0.0
    gamma = 0.000001832 ± 0.000011204,
#    Pbdot = -0.0975723e-12 - 0.0018e-12 ± 0.0463125e-12,
    Pbdot = -0.0975723e-12 ± 0.0463125e-12,
    s = 0.992298 ± 0.023784,
    m2_shapiro = 0.384013 ± 0.346765,
    m1_init = 1.241308,
    m2_init = 0.950913
    )    

obs_params_dataset["J1952+2630_DDFWHE"] = ObsParams(
    Pb = 0.39188061942246 ± 0.00000030769805,
    T0 = 55407.6182211277 ± 0.0022697669,
    e0 = 0.0000409895 ± 0.0000013922,
    omega0 = 289.961323343612 ± 2.085070587999,
    x0 = 2.798183211 ± 0.000007394,
    omegadot = 1.7009239 ± 0.2634653,
#   gamma = GAMMA          0.000000225 0.0
    gamma = 0.000001824 ± 0.000011215,
#    Pbdot = -0.1007306e-12 + 0.0018000e-12 ± 0.0463139e-12,
    Pbdot = -0.1007306e-12 ± 0.0463139e-12,
    h3 = 0.000001305397 ± 0.000000591883,
#    varsigma = 0.883989 ± 0.175774,
    m1_init = 1.241308,
    m2_init = 0.950913
    )

obs_params_dataset["J1952+2630_DDSTG"] = ObsParams(
    Pb = 0.39188049478442 ± 0.00000021533086,
    T0 = 55407.6178241798 ± 0.0017828332,
    e0 = 0.0000422372 ± 0.0000010009,
    omega0 = 289.596685703459 ± 1.637720152256,
    x0 = 2.798177056 ± 0.000001905,
    omegadot = 1.5942128 ± 0.2634732,
#   gamma = GAMMA          0.000000225 0.0
#    gamma = 0.000001824 ± 0.000011215,
#    Pbdot = -0.0971310e-12 - 0.0018e-12 ± 0.0463139e-12,
    Pbdot = -0.0935800e-12 ± 0.001e-12,
    h3 = 0.00000173760 ± 0.0000001,
    varsigma = 0.7192647 ± 0.01,
    m1_init = 1.241308,
    m2_init = 0.950913
    )

obs_params_dataset["J1952+2630_DDFWHE_23-26_efac"] = ObsParams(
    Pb = 0.39188054006575 ± 0.00000013384599,
    T0 = 55407.6164803160 ± 0.0017013006,
    e0 = 0.0000430302 ± 0.0000008502,
    omega0 = 288.361953619811 ± 1.562855646502,
    x0 = 2.798167234 ± 0.000008294,
    omegadot = 1.6329910 ± 0.1145993,
    gamma = 0.000006818 ± 0.000005028,
#    Pbdot = -0.0933166e-12 - 0.0018e-12 ± 0.0032936e-12,
    Pbdot = -0.0933166e-12 ± 0.0032936e-12,
    h3 = 0.000001934141 ± 0.000000248389,
    varsigma = 0.645961 ± 0.105603,
    m1_init = 1.241308,
    m2_init = 0.950913
    )

obs_params_dataset["J1952+2630_DDFWHE_23-39_efac"] = ObsParams(
    Pb = 0.39188049109643 ± 0.00000005167042,
    T0 = 55407.6186090237 ± 0.0009775777,
    e0 = 0.0000419971 ± 0.0000004048,
    omega0 = 290.317731580282 ± 0.898029728255,
    x0 = 2.798181164 ± 0.000003802,
    omegadot = 1.5910538 ± 0.0442405,
#    gamma = -0.000002244 ± 0.000002389,
#    Pbdot = -0.0936739e-12 - 0.0018e-12 ± 0.0006421e-12,
    Pbdot = -0.0936739e-12 ± 0.0006421e-12,
    h3 = 0.000001533496 ± 0.000000179930,
    varsigma = 0.734309 ± 0.079082,
    m1_init = 1.241308,
    m2_init = 0.950913
    )

obs_params_dataset["J1952+2630_DDFWHE_4o_32"] = ObsParams(
    Pb = 0.39188051335950 ± 0.00000004222182,
    T0 = 55407.6172166659 ± 0.0006879391,
    e0 = 0.0000425971 ± 0.0000002225,
    omega0 = 289.038525212031 ± 0.631959273535,
    x0 = 2.798173891 ± 0.000002358,
    omegadot = 1.6101278 ± 0.0361506,
#    gamma = 0.000000654 ± 0.000001884,
#    Pbdot = -0.0950446e-12 + 0.0018e-12 ± 0.0007707e-12,
    Pbdot = -0.0937366e-12 ± 0.0007707e-12,
    h3 = 0.000001952196 ± 0.000000082282,
    varsigma = 0.690174 ± 0.031625,
    m1_init = 1.241308,
    m2_init = 0.950913
    )

obs_params_dataset["J1952+2630_DDFWHE_4o_32_1SXPBDOT"] = ObsParams(
    Pb = 0.39188051335950 ± 0.00000004222182,
    T0 = 55407.6172166659 ± 0.0006879391,
    e0 = 0.0000425971 ± 0.0000002225,
    omega0 = 289.038525212031 ± 0.631959273535,
    x0 = 2.798173891 ± 0.000002358,
    omegadot = 1.6101278 ± 0.0361506,
#    gamma = 0.000000654 ± 0.000001884,
#    Pbdot = -0.0950446e-12 + 0.0018e-12 ± 0.0007707e-12,
    Pbdot = -0.0937366e-12 - 0.005e-12 ± 0.0007707e-12,
    h3 = 0.000001952196 ± 0.000000082282,
    varsigma = 0.690174 ± 0.031625,
    m1_init = 1.241308,
    m2_init = 0.950913
    )

obs_params_dataset["J1952+2630_DDFWHE_4o_32_1.5M_1SXPBDOT"] = ObsParams(
    Pb = 0.39188076252287 ± 0.00000004596516,
    T0 = 55407.6163956321 ± 0.0007504388,
    e0 = 0.0000422904 ± 0.0000003866,
    omega0 = 288.284436112367 ± 0.689374033012,
    x0 = 2.798170152 ± 0.000006218,
    omegadot = 1.8233860 ± 0.0393557,
#    gamma = 0.000001247 ± 0.000001747,
#    Pbdot = -0.1350397e-12 - 0.0018e-12 ± 0.0007718e-12,
    Pbdot = -0.1350397e-12 - 0.0051e-12 ± 0.0007718e-12,
    h3 = 0.000001179195 ± 0.000000081866,
    varsigma = 0.483049 ± 0.066909,
    m1_init = 1.241308,
    m2_init = 0.950913
    )

    #=
obs_params_dataset["J1952+2630_DDFWHE_4o_32"] = ObsParams(
    Pb = 0.39188051738605 ± 0.00000004240271,
    T0 = 55407.6169782060 ± 0.0006901094,
    e0 = 0.0000423401 ± 0.0000002202,
    omega0 = 288.819546585191 ± 0.633952970133,
    x0 = 2.798175542 ± 0.000002313,
    omegadot = 1.6135625 ± 0.0363054,
#    gamma = 0.000000269 ± 0.000001875,
#    Pbdot = -0.0939152e-12 - 0.0018e-12 ± 0.0007707e-12,
    Pbdot = -0.0939152e-12 ± 0.0007707e-12,
    h3 = 0.000001811314 ± 0.000000082192,
    varsigma = 0.694783 ± 0.033725,
    m1_init = 1.241308,
    m2_init = 0.950913
    )
=#

#=
obs_params_dataset["J1952+2630_DDFWHE_4o_32_XPBDOT_1S"] = ObsParams(
    Pb = 0.39188051738605 ± 0.00000004240271,
    T0 = 55407.6169782060 ± 0.0006901094,
    e0 = 0.0000423401 ± 0.0000002202,
    omega0 = 288.819546585191 ± 0.633952970133,
    x0 = 2.798175542 ± 0.000002313,
    omegadot = 1.6135625 ± 0.0363054,
#    gamma = 0.000000269 ± 0.000001875,
#    Pbdot = -0.0939152e-12 - 0.0018e-12 ± 0.0007707e-12,
    Pbdot = -0.0939152e-12 - 0.005e-12 ± 0.0007707e-12,
    h3 = 0.000001811314 ± 0.000000082192,
    varsigma = 0.694783 ± 0.033725,
    m1_init = 1.241308,
    m2_init = 0.950913
    )
=#

obs_params_dataset["J1952+2630_DDFWHE_4o_32_XPBDOT"] = ObsParams(
    Pb = 0.39188051738605 ± 0.00000004240271,
    T0 = 55407.6169782060 ± 0.0006901094,
    e0 = 0.0000423401 ± 0.0000002202,
    omega0 = 288.819546585191 ± 0.633952970133,
    x0 = 2.798175542 ± 0.000002313,
    omegadot = 1.6135625 ± 0.0363054,
#    gamma = 0.000000269 ± 0.000001875,
#    Pbdot = -0.0939152e-12 - 0.0018e-12 ± 0.0007707e-12,
    Pbdot = -0.0939152e-12 ± 0.005e-12,
    h3 = 0.000001811314 ± 0.000000082192,
    varsigma = 0.694783  ± 0.033725,
    m1_init = 1.241308,
    m2_init = 0.950913
    )


#=    
obs_params_dataset["J1952+2630_DDFWHE_23-39"] = ObsParams(
    Pb = 0.39188047365222 ± 0.00000005079697,
    T0 = 55407.6184824073 ± 0.0009435680,
    e0 = 0.0000423088 ± 0.0000003998,
    omega0 = 290.201314405219 ± 0.866787508168,
    x0 = 2.798178093 ± 0.000003454,
    omegadot = 1.5761222 ± 0.0434927,
    gamma = 0.000000637 ± 0.000002355,
#    Pbdot = -0.0933530e-12 - 0.0018e-12 ± 0.0006284e-12,
    Pbdot = -0.0933530e-12 ± 0.0006284e-12,
    h3 = 0.000001785735 ± 0.000000186858,
    varsigma = 0.766822 ± 0.064952,
    m1_init = 1.241308,
    m2_init = 0.950913
    )
=#

#=   
obs_params_dataset["J2222−0137_Guo"] = ObsParams(
    Pb = 2.445759995471 ± 0.000000000006,
    T0= 58001.200912600 ± 0.000000002,
    e0 = 0.0003809326 ± 0.0000000138,
    omega0 = 120.45820 ± 0.00143,
    x0 = 10.8480213 ± 0.0000002,
    omegadot = 0.09607 ± 0.00048,
#   Pbdot = -0.0143e-12 ± 0.0076e-12,
    Pbdot = -0.0080267e-12 ± 0.0076e-12,
    m2_shapiro = 1.3153 ± 0.009,
    s = 0.996623 ± 0.000115
    )
=#

obs_params_dataset["J2222−0137_Guo_ELL1H+"] = ObsParams(
    Pb = 2.445759995471 ± 0.000000000006,
    T0 = 58001.200912600 ± 0.000000002,
    eps1 = 0.00032836 ± 0.00000002,
    eps2 = -0.000193092 ± 0.000000009,
    x0 = 10.8480213 ± 0.0000002,
    omegadot = 0.09607 ± 0.00048,
#    Pbdot = 0.2554e-12 ± 0.0074e-12,
    Pbdot = -0.0098e-12 ± 0.0078e-12,
    h3 = 5.052e-6 ± 0.027e-6,
    varsigma = 0.9212 ± 0.0014,
    m1_init = 1.81,
    m2_init = 1.31
    )

obs_params_dataset["J2222−0137_Guo_DDK"] = ObsParams(
    Pb = 2.44576437 ± 0.00000002,
    T0 = 58002.01928 ± 0.00001,
    e0 = 0.00038092 ± 0.00000001,
    omega0 = 120.458 ± 0.002,
    x0 = 10.8480235 ± 0.0000002,
    omegadot = 0.09605 ± 0.00048,
#    Pbdot = 0.2509e-12 ± 0.0074e-12,
    Pbdot = -0.0143e-12 ± 0.0078e-12,
    m2_shapiro = 1.315 ± 0.012,
    s = 0.99661 ± 0.00012,
    m1_init = 1.81,
    m2_init = 1.31
    )

obs_params_dataset["J2222−0137_DDK"] = ObsParams(
    Pb = 2.44576437112767 ± 0.00000002221810,
    T0 = 58002.0192842713 ± 0.0000112869,
    e0 = 0.0003809111 ± 0.0000000133,
    omega0 = 120.459010621444 ± 0.001661336315,
    x0 = 10.848023346 ± 0.000000176,
    omegadot = 0.0961852 ± 0.0004883,
#    Pbdot = 0.2509e-12 ± 0.0074e-12,
    Pbdot = -0.0115657e-12 ± 0.0081135e-12,
    m2_shapiro = 1.329476 ± 0.010868,
    s = 0.99651 ± 0.00012,
    m1_init = 1.81,
    m2_init = 1.31
    )

obs_params_dataset["J2222−0137_DDK_XDOT"] = ObsParams(
    Pb = 2.44576437993512 ± 0.00000002233131,
    T0 = 58002.0192752407 ± 0.0000114774,
    e0 = 0.0003809287 ± 0.0000000140,
    omega0 = 120.457681353810 ± 0.001689384445,
    x0 = 10.848023632 ± 0.000000188,
    omegadot = 0.0963789 ± 0.0004908,
#    Pbdot = 0.2509e-12 ± 0.0074e-12,
    Pbdot = -0.0161612e-12 ± 0.0081452e-12,
    m2_shapiro = 1.309790 ± 0.011772,
    s = 0.99669 ± 0.00013,
    m1_init = 1.81,
    m2_init = 1.31
    )

obs_params_dataset["J2222−0137_DD"] = ObsParams(
    Pb = 2.44576438024568 ± 0.00000002232699,
    T0 = 58002.0192754945 ± 0.0000114642,
    e0 = 0.0003809317 ± 0.0000000138,
    omega0 = 120.457718993631 ± 0.001687436748,
    x0 = 10.848023666 ± 0.000000187,
    omegadot = 0.0963854 ± 0.0004908,
#    Pbdot = -0.0120165e-12 + 0.2652000e-12 ± 0.0074808e-12,
    Pbdot = -0.0120183e-12 ± 0.0074808e-12,
#    Pbdot = -0.0080654e-12 ± 0.0074808e-12,
    m2_shapiro = 1.307742 ± 0.011662,
    s = 0.996706 ± 0.000125,
    m1_init = 1.8223,
    m2_init = 1.3162
    )    


obs_params_dataset["J2222−0137_DD_bestXPBDOT"] = ObsParams(
    Pb = 2.44576438024778 ± 0.00000002232699,
    T0 = 58002.0192754937 ± 0.0000114642,
    e0 = 0.0003809317 ± 0.0000000138,
    omega0 = 120.457718872219 ± 0.001687437091,
    x0 = 10.848023666 ± 0.000000187,
    omegadot = 0.0963855 ± 0.0004908,
#    Pbdot = -0.0080019e-12 + 0.2611851e-12 ± 0.0074808e-12,
    Pbdot = -0.0080019e-12 ± 0.0074808e-12,
#    Pbdot = -0.0080654e-12 ± 0.0074808e-12,
    m2_shapiro = 1.307741 ± 0.011662,
    s = 0.996706 ± 0.000125,
    m1_init = 1.8223,
    m2_init = 1.3162
    )    

obs_params_dataset["J2222−0137_DDFWHE"] = ObsParams(
    Pb = 2.44576438025421 ± 0.00000002232698,
    T0 = 58002.0192754929 ± 0.0000114642,
    e0 = 0.0003809317 ± 0.0000000138,
    omega0 = 120.457718692649 ± 0.001687436968,
    x0 = 10.848023666 ± 0.000000187,
    omegadot = 0.0963855 ± 0.0004908,
#    Pbdot = -0.0120165e-12 + 0.2652000e-12 ± 0.0074808e-12,
    Pbdot = -0.0120151e-12 ± 0.0074808e-12,
#    Pbdot = -0.0080654e-12 ± 0.0074808e-12,
    h3 = 0.000005047484 ± 0.000000027016,
    varsigma = 0.921937 ± 0.001423,
    m1_init = 1.8223,
    m2_init = 1.3162
    )

obs_params_dataset["J2222−0137_DDFWHE_bestXPBDOT"] = ObsParams(
    Pb = 2.44576438025448 ± 0.00000002232699,
    T0 = 58002.0192754933 ± 0.0000114642,
    e0 = 0.0003809317 ± 0.0000000138,
    omega0 = 120.457718811164 ± 0.001687436998,
    x0 = 10.848023666 ± 0.000000187,
    omegadot = 0.0963856 ± 0.0004908,
#    Pbdot = -0.0120165e-12 + 0.2652000e-12 ± 0.0074808e-12,
    Pbdot = -0.0080009e-12 ± 0.0074808e-12,
#    Pbdot = -0.0080654e-12 ± 0.0074808e-12,
    h3 = 0.000005047482 ± 0.000000027016,
    varsigma = 0.921937 ± 0.001423,
    m1_init = 1.8223,
    m2_init = 1.3162
    )

obs_params_dataset["J2222−0137_sim_DD"] = ObsParams(
    Pb = 2.44576436212514 ± 0.00000000144788,
    T0 = 58002.0192821774 ± 0.0000019632,
    e0 = 0.0003809227 ± 0.0000000011,
    omega0 = 120.458702637165 ± 0.000288970653,
    x0 = 10.848023540 ± 0.000000017,
    omegadot = 0.0959873 ± 0.0000318,
#    Pbdot = -0.0080018e-12 + 0.2633700e-12 ± 0.0004217e-12,
    Pbdot = -0.0080018e-12 ± 0.0004217e-12,
    m2_shapiro = 1.315034 ± 0.000945,
    s = 0.996624 ± 0.000011,
    m1_init = 1.8223,
    m2_init = 1.3162
    )    

obs_params_dataset["J1738+0333"] = ObsParams(
    Pb = 0.3547907398724 ± 0.0000000000013,
    e0 = 0.00000034 ± 0.00000011,
    x0 = 0.343429130 ± 0.0000017,
#    Pbdot = 0.2509e-12 ± 0.0074e-12,
    Pbdot = -0.0259e-12 ± 0.0033e-12,
    m2 = 0.181 ± 0.008,
    q = 8.1 ± 0.2
    )

obs_params_dataset["J1738+0333_Guo"] = ObsParams(
    Pb = 0.3547907398724 ± 0.0000000000013,
    e0 = 0.00000034 ± 0.00000011,
#    Pbdot = 0.2509e-12 ± 0.0074e-12,
    Pbdot = -0.02772e-12 ± -0.00064e-12,
    m2 = 0.1817 ± 0.0073,
    q = 8.1 ± 0.2
    )


obs_params_dataset["Triple System"] = ObsParams(
    deltaN = 0.5 ± 0.9
    )

obs_params_dataset["Double Pulsar"] = ObsParams(
    Pb = 0.1022515592973 ± 0.0000000000010,
    T0 = 55700.233017540 ± 0.000000013,
    e0 = 0.087777023 ± 0.000000061,
    omega0 = 204.753686 ± 0.000047,
    x0 = 1.415028603 ± 0.000000092,
    omegadot = 16.899323 ± 0.000013,
    gamma = 0.384045e-3 ± 0.000094e-3,
    Pbdot = -1.247752e-12 ± 0.000079e-12,
    r = 6.162e-6 ± 0.021e-6,
    s = 0.999936 ± 0.000010,
    m1_init = 1.338185,
    m2_init = 1.248868
    )
