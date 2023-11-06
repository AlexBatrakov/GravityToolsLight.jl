#--------------------------------------------------------------------------------------------------------------
# Общие настройки Tempo, которые могут включать различные расширенные настройки
mutable struct GeneralTempoSettings{T <: AbstractTempoVersion}
    basic_settings::BasicTempoSettings{T}
    global_iter_settings::Union{GlobalIterationsSettings, Nothing}
    param_sweep_settings::Union{ParameterSweepSettings, Nothing}
end

# Конструкторы
function GeneralTempoSettings(
    basic_settings::BasicTempoSettings{T}, 
    global_iter_settings::Union{GlobalIterationsSettings, Nothing} = nothing,
    param_sweep_settings::Union{ParameterSweepSettings, Nothing} = nothing
    ) where {T <: AbstractTempoVersion}
    return GeneralTempoSettings{T}(basic_settings, global_iter_settings, param_sweep_settings)
end

#--------------------------------------------------------------------------------------------------------------

mutable struct GeneralTempoFramework{T1 <: TestParameters, T2 <: RefinementSettings}
    tsets::GeneralTempoSettings
    test_params::T1
    ref_sets::T2
    grid::AdaptiveRefinement2DGrid
end

function Base.show(io::IO, tf::GeneralTempoFramework)
    println(io, "Tempo framework:")
    println(io, tf.tsets)
    println(io, tf.test_params)
    print(io, tf.ref_sets)
	return nothing
end


# function GeneralTempoFramework(test::GeneralTest, obs_params::ObsParams, gsets::GridSetttings)
#     param1_grid = collect(LinRange(test.param1.min, test.param1.max, test.param1.N))
#     param2_grid = collect(LinRange(test.param2.min, test.param2.max, test.param2.N))
#     grid = Refinement2DGrid(Dict(), param1_grid, param2_grid)
#     return PKFramework(test, obs_params, gsets, grid)
# end

GeneralTempoFramework(;tsets::GeneralTempoSettings, test_params::T1, ref_sets::T2) where {T1 <: TestParameters, T2 <: RefinementSettings} = TempoFramework{T1, T2}(tsets, test_params, ref_sets)

function GeneralTempoFramework(tsets::GeneralTempoSettings, test_params::TestParameters, ref_sets::RefinementSettings)
    grid = AdaptiveRefinement2DGrid(test_params.x, test_params.y, ref_sets)
    return GeneralTempoFramework(tsets, test_params, ref_sets, grid)
end