KObsType = NamedTuple{K_list, NTuple{length(K_list), Measurement{Float64}}}
PKObsType = NamedTuple{PK_list, NTuple{length(PK_list), Measurement{Float64}}}
XObsType = NamedTuple{X_list, NTuple{length(X_list), Measurement{Float64}}}

mutable struct ClassicalTest
    K_params_obs::KObsType
    PK_params_obs::PKObsType
    X_params_obs::XObsType
    grid::SimpleGrid
    N_refinement::Int64
    CL::Float64
    function ClassicalTest(K_params_obs, PK_params_obs)
        return new(K_params_obs, PK_params_obs, XObsType(zeros(length(X_list))))
    end
    function ClassicalTest(K_params_obs, PK_params_obs, X_params_obs)
        return new(K_params_obs, PK_params_obs, X_params_obs)
    end
end

function Base.show(io::IO, cl::ClassicalTest)
    println(io, "ClassicalTest")
	println(io, "Observed Keplerian parameters: ", cl.K_params_obs)
	println(io, "Observed Post-Keplerian parameters: ", cl.PK_params_obs)
	println(io, "Observed Extra parameters: ", cl.X_params_obs)
	println(io, "Desired refinement level: ", cl.N_refinement)
    println(io, "Desired confidence level: ", cl.CL)
	return nothing
end

function find_initial_masses(cl::ClassicalTest, pf::DEFPhysicalFramework)
    PK_first = :k
    PK_second = :gamma
    for PK in keys(cl.PK_params_obs)
        if abs(cl.PK_params_obs[PK].val / cl.PK_params_obs[PK].err) > abs(cl.PK_params_obs[PK_first].val / cl.PK_params_obs[PK_first].err)
            PK_first = PK
            PK_second = PK_first
        elseif abs(cl.PK_params_obs[PK].val / cl.PK_params_obs[PK].err) > abs(cl.PK_params_obs[PK_second].val / cl.PK_params_obs[PK_second].err) && PK != PK_first
            PK_second = PK
        end
    end

    println(PK_first, " ", PK_second)

    function find_intersection!(F, x)
        pf.bnsys.psr.mass = x[1]
        pf.bnsys.comp.mass = x[2]
        interpolate_bnsys!(pf)
        F[1] = (pf.bnsys.PK_params[PK_first] / cl.PK_params_obs[PK_first].val) - 1.0 #/ cl.PK_params_sigma[PK_first]
        F[2] = (pf.bnsys.PK_params[PK_second] / cl.PK_params_obs[PK_second].val) - 1.0 #/ cl.PK_params_sigma[PK_first]
    end
    solution = nlsolve(find_intersection!, [pf.bnsys.psr.mass; pf.bnsys.comp.mass])
#    println(solution)
    m1, m2 = solution.zero
    pf.bnsys.psr.mass = m1
    pf.bnsys.comp.mass = m2
    return (m1 = m1, m2 = m2)
end

function get_chisqr(cl::ClassicalTest, pf::DEFPhysicalFramework)
    chisqr = 0.0
    for PK in keys(cl.PK_params_obs)
        if cl.PK_params_obs[PK].err != 0.0
            chisqr += ((pf.bnsys.PK_params[PK] - cl.PK_params_obs[PK].val) / cl.PK_params_obs[PK].err)^2
        end
    end
    for X in keys(cl.X_params_obs)
        if cl.X_params_obs[X].err != 0.0
            chisqr += ((pf.bnsys.X_params[X] - cl.X_params_obs[X].val) / cl.X_params_obs[X].err)^2
        end
    end
    return chisqr
end

function check_terms_in_chisqr(cl::ClassicalTest, pf::DEFPhysicalFramework)
    println("(theory - observation) / sigma")
    for PK in keys(cl.PK_params_obs)
        if cl.PK_params_obs[PK].err != 0.0
            println("$PK: ", ((pf.bnsys.PK_params[PK] - cl.PK_params_obs[PK].val) / cl.PK_params_obs[PK].err))
        end
    end
    for X in keys(cl.X_params_obs)
        if cl.X_params_obs[X].err != 0.0
            println("$X: ", ((pf.bnsys.X_params[X] - cl.X_params_obs[X].val) / cl.X_params_obs[X].err))
        end
    end
end

function find_best_masses(cl::ClassicalTest, pf::DEFPhysicalFramework)
    function get_chisqr_local(x)
        pf.bnsys.psr.mass = x[1]
        pf.bnsys.comp.mass = x[2]
        interpolate_bnsys!(pf)
        return get_chisqr(cl::ClassicalTest, pf::DEFPhysicalFramework)
    end
    sol = optimize(get_chisqr_local, [pf.bnsys.psr.mass, pf.bnsys.comp.mass])
    m1, m2 = Optim.minimizer(sol)
    pf.bnsys.psr.mass, pf.bnsys.comp.mass = m1, m2
    return (m1 = m1, m2 = m2)
end


function calculate!(cl::ClassicalTest, pf::DEFPhysicalFramework)

    function calculate_chisqr(log10alpha0::Float64, beta0::Float64)



        pf.bnsys.psr.mass = 1.5
        pf.bnsys.comp.mass = 1.3
        interpolate_bnsys!(pf)
        return (chisqr = chisqr)
    end

    lvl = quantile(Chisq(2), cl.CL)

    function contour_cell_selector(i_cell::Int64, j_cell::Int64, grid::SimpleGrid)
        chisqr_min = grid.params[:chisqr_min]
        chisqr_cell = @view grid.value[:chisqr][i_cell:i_cell+1,j_cell:j_cell+1]
        chisqr_cell_min = minimum(chisqr_cell)
        chisqr_cell_max = maximum(chisqr_cell)
        return chisqr_cell_min < chisqr_min + lvl < chisqr_cell_max
    end

    function calculate_params!(grid::SimpleGrid)
        grid.params[:chisqr_min] = min(minimum(grid.value[:chisqr]), grid.params[:chisqr_min])
        return nothing
    end

    precalculate_Grid(grid, calculate_chisqr, calc_params!)
    for i in 1:cl.N_refinement
        cl.grid = refine_Grid(cl.grid, calculate_chisqr, contour_cell_selector, calculate_params!)
    end
end