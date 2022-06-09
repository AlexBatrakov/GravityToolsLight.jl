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
    CL::Float64
    plot_type::String
end

function Base.show(io::IO, gsets::GridSetttings)
    println(io, "Grid settings:")
	println(io, "   Desired refinement level: ", gsets.N_refinement)
    println(io, "   Desired confidence level: ", gsets.CL)
    print(io,   "   Type of the plot: ", gsets.plot_type)
	return nothing
end

GridSetttings(;N_refinement, CL, plot_type) = GridSetttings(N_refinement, CL, plot_type)

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