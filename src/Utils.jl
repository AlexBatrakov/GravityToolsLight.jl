lvl_1σ = quantile(Chisq(2), 0.682689492137086)
lvl_2σ = quantile(Chisq(2), 0.954499736103642)
lvl_3σ = quantile(Chisq(2), 0.997300203936740)
lvl_4σ = quantile(Chisq(2), 0.999936657516334)
lvl_5σ = quantile(Chisq(2), 0.999999426696856)
lvl_6σ = quantile(Chisq(2), 0.999999998026825)
lvl_7σ = quantile(Chisq(2), 0.999999999997440)

lvl_68CL = quantile(Chisq(2), 0.68)
lvl_90CL = quantile(Chisq(2), 0.90)
lvl_95CL = quantile(Chisq(2), 0.95)
lvl_99CL = quantile(Chisq(2), 0.99)
lvl_997CL = quantile(Chisq(2), 0.997)

abstract type AbstractGravityTest end

Range = NamedTuple{(:min, :max, :N), Tuple{Float64, Float64, Int64}}
RangedParameter = NamedTuple{(:name, :min, :max, :N), Tuple{String, Float64, Float64, Int64}}


struct GeneralTest{T1 <: Union{Float64, Range}, T2 <: Union{Float64, Range}} <: AbstractGravityTest
    psrname::String
    eosname::String
    alpha0::T1
    log10alpha0::T1
    beta0::T2
    param1::RangedParameter
    param2::RangedParameter
end

GeneralTest(test::GeneralTest, eosname) = GeneralTest(test.psrname, eosname, test.alpha0, test.log10alpha0, test.beta0, test.param1, test.param2) 

function Base.show(io::IO, test::GeneralTest)
    println(io, "General test:")
    println(io,     "   Name of the pulsar:  ", test.psrname)
    if test.eosname != ""
        println(io, "   Equation of state:   ", test.eosname)
    end
    if typeof(test.alpha0) == Float64 && !isnan(test.alpha0)
        println(io, "   alpha0:              ", test.alpha0)
        println(io, "   log10|alpha0|:       ", test.log10alpha0)
    elseif typeof(test.alpha0) == Range
        println(io, "   alpha0 range:        ", test.alpha0)
        println(io, "   log10|alpha0| range: ", test.log10alpha0)
    end
    if typeof(test.beta0) == Float64 && !isnan(test.beta0)
        println(io, "   beta0:               ", test.beta0)
    elseif typeof(test.beta0) == Range
        println(io, "   beta0 range:         ", test.beta0)
    end
    println(io,     "   First parameter:     ", test.param1)
    print(io,       "   Second parameter:    ", test.param2)
	return nothing
end

function read_param!(param, alpha0, log10alpha0, beta0)
    if param.name == "alpha0"
        alpha0 = (min = param.min, max = param.max, N = param.N)
        log10alpha0 = (min = log10(abs(param.min)), max = log10(abs(param.max)), N = param.N)
    elseif param.name == "log10alpha0"
        log10alpha0 = (min = param.min, max = param.max, N = param.N)
        alpha0 = (min = -exp10(param.min), max = -exp10(param.max), N = param.N)
    elseif param.name == "beta0"
        beta0 = (min = param.min, max = param.max, N = param.N)
    end
    return alpha0, log10alpha0, beta0
end

function GeneralTest(;psrname::String, eosname::String="", alpha0=NaN, log10alpha0=NaN, beta0=NaN, param1, param2)
    if !isnan(alpha0)
        log10alpha0 = log10(abs(alpha0))
    elseif !isnan(log10alpha0)
        log10alpha0 = Float64(log10alpha0)
        alpha0 = -exp10(log10alpha0)
    end
    alpha0, log10alpha0, beta0 = read_param!(param1, alpha0, log10alpha0, beta0)
    alpha0, log10alpha0, beta0 = read_param!(param2, alpha0, log10alpha0, beta0)
    return GeneralTest(psrname, eosname, alpha0, log10alpha0, beta0, param1, param2)
end

struct GridSetttings
    N_refinement::Int64
    CL::Vector{Float64}
    contours::Vector{Float64}
    refinement_type::String
    delta_chisqr_max::Float64
    delta_chisqr_diff::Float64
    gr_in_chisqr::Bool
    function GridSetttings(N_refinement, CL, contours, refinement_type, delta_chisqr_max, delta_chisqr_diff, gr_in_chisqr)
        CL = typeof(CL) == Float64 ? [CL] : CL
        return new(N_refinement, isempty(CL) ? cdf(Chisq(2), contours) : CL, isempty(contours) ? quantile.(Chisq(2), CL) : contours, refinement_type, delta_chisqr_max, delta_chisqr_diff, gr_in_chisqr)
    end
end

function Base.show(io::IO, gsets::GridSetttings)
    println(io, "Grid settings:")
	println(io, "   Desired refinement level: ", gsets.N_refinement)
    println(io, "   Desired confidence levels: ", gsets.CL)
    println(io, "   Desired Δχ2 contours: ", gsets.contours)
    println(io, "   Type of the refinements: ", gsets.refinement_type)
    println(io, "   Maximum Δχ2 value: ", gsets.delta_chisqr_max)
    println(io, "   Maximum difference in Δχ2 value: ", gsets.delta_chisqr_diff)
    print(io,   "   Include GR in Δχ2 value: ", gsets.gr_in_chisqr)
	return nothing
end


GridSetttings(;N_refinement=1, CL=Float64[], contours=Float64[], refinement_type="nice", delta_chisqr_max=10.0, delta_chisqr_diff=1.0, gr_in_chisqr=false) = GridSetttings(N_refinement, CL, contours, refinement_type, delta_chisqr_max, delta_chisqr_diff, gr_in_chisqr)


#=
function get_label(name)
    if name == "PBDOT"
        return L"\dot{P}_\mathrm{b}\, (10^{-12})"
    elseif name == "M2"
        return L"m_{\mathrm{c}}\,(M_\odot)"
    elseif name == "MTOT"
        return L"m_{\mathrm{tot}}\,(M_\odot)"
    elseif name == "GAMMA"
        return L"\gamma"
    elseif name == "XDOT"
        return L"\dot{x}\,(10^{-12} s/s)"
    elseif name == "OMDOT"
        return L"\dot{\omega}"
    elseif name == "COSI"
        return L"\cos i"
    elseif name == "DTHETA"
        return L"\delta_\theta"
    elseif name == "log10alpha0"
        return L"\log_{10}|\alpha_0|"
    elseif name == "beta0"
        return L"\beta_0"
    elseif name == "H3"
        return L"h_3"
    elseif name == "VARSIGMA"
        return L"\varsigma"
    end
end
=#