mutable struct ClassicalTest
    K_params::KType
    PK_params_obs::PKType
    PK_params_sigma::PKType
#    PK_params_min::PKType
    grid::SimpleGrid
    N_refinement::Int64
    CL::Float64
    function ClassicalTest(K_params::KType, PK_params_obs::PKType, PK_params_sigma::PKType)
        return new(K_params, PK_params_obs, PK_params_sigma)
    end
end

function find_initial_masses(cl::ClassicalTest, pf::DEFPhysicalFramework)
    PK_first = :k
    PK_second = :gamma
    for PK in keys(cl.PK_params_obs)
        if abs(cl.PK_params_obs[PK] / cl.PK_params_sigma[PK]) > abs(cl.PK_params_obs[PK_first] / cl.PK_params_sigma[PK_first])
            PK_first = PK
            PK_second = PK_first
        elseif abs(cl.PK_params_obs[PK] / cl.PK_params_sigma[PK]) > abs(cl.PK_params_obs[PK_second] / cl.PK_params_sigma[PK_second]) && PK != PK_first
            PK_second = PK
        end
    end

    println(PK_first, " ", PK_second)

    function find_intersection!(F, x)
        pf.bnsys.psr.mass = x[1]
        pf.bnsys.comp.mass = x[2]
        interpolate_bnsys!(pf)
        F[1] = (pf.bnsys.PK_params[PK_first] / cl.PK_params_obs[PK_first]) - 1.0 #/ cl.PK_params_sigma[PK_first]
        F[2] = (pf.bnsys.PK_params[PK_second] / cl.PK_params_obs[PK_second]) - 1.0 #/ cl.PK_params_sigma[PK_first]
    end
    solution = nlsolve(find_intersection!, [pf.bnsys.psr.mass; pf.bnsys.comp.mass])
    println(solution)
    m1, m2 = solution.zero
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