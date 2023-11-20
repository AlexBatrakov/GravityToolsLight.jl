#--------------------------------------------------------------------------------------------------------------

"""
    struct BasicTempoKeys

Holds the key settings for running the Tempo software.

# Fields
- `silent::Bool`: If true, Tempo runs without printing progress to the console.
- `print_output::Bool`: If true, the output of Tempo is printed.
- `iterative_mode::Bool`: If true, Tempo runs in iterative mode.
- `fit_EFACs_EQUADs::Bool`: If true, EFAC and EQUAD parameters are fitted.

# Constructor
    BasicTempoKeys(;silent = true, print_output = false, iterative_mode = true, fit_EFACs_EQUADs = false)
"""
struct BasicTempoKeys
    silent::Bool
    print_output::Bool
    save_internal_iterations::Bool
    iterative_mode::Bool
    fit_EFACs_EQUADs::Bool
    BasicTempoKeys(;silent = true, print_output = false, save_internal_iterations=false, iterative_mode=true, fit_EFACs_EQUADs::Bool = false) = new(silent, print_output, save_internal_iterations, iterative_mode, fit_EFACs_EQUADs)
end

# Implementation of the show method for BasicTempoKeys
function Base.show(io::IO, keys::BasicTempoKeys)
    println(io, "Tempo keys:")
    print(io, "        Silent mode: ", keys.fit_EFACs_EQUADs)
    print(io, "        Fit EFAC and EQUAD parameters: ", keys.fit_EFACs_EQUADs)
    return nothing
end

#--------------------------------------------------------------------------------------------------------------

"""
    struct BasicTempoSettings{T <: AbstractTempoVersion}

Holds all the basic settings required to run a single instance of the Tempo software.

# Fields
- `work_dir::String`: The working directory for the Tempo run.
- `version::T`: The version of Tempo to be used, must be a subtype of `AbstractTempoVersion`.
- `par_file_init::String`: The initial parameter file for Tempo.
- `tim_file::String`: The timing data file for Tempo.
- `flags::String`: Additional command line flags for Tempo.
- `keys::BasicTempoKeys`: Key settings for running Tempo.
- `tparams::Vector{GeneralTempoParameter}`: A vector of tempo parameters to be used.

# Constructor
    BasicTempoSettings(;work_dir, version, par_file_init, tim_file, flags = "", keys = BasicTempoKeys(), tparams = GeneralTempoParameter[])
"""
struct BasicTempoSettings{T <: AbstractTempoVersion}
    work_dir::String
    version::T
    par_file_init::String
    tim_file::String
    flags::String
    keys::BasicTempoKeys
    tparams::Vector{GeneralTempoParameter}
end

# Implementation of the show method for BasicTempoSettings
function Base.show(io::IO, bsets::BasicTempoSettings)
    println(io, "Basic Tempo settings:")
    println(io, "   Working directory: ", bsets.work_dir)
    println(io, "   Version: ", bsets.version)
    println(io, "   Initial par file: ", bsets.par_file_init)
    println(io, "   Working tim file: ", bsets.tim_file)
    println(io, "   Selected additional flags: ", bsets.flags)
    println(io, "   ", bsets.keys)
    println(io, "   Tempo parameters used:", bsets.tparams)
    return nothing
end

# Constructor implementation
BasicTempoSettings(;
    work_dir,
    version,
    par_file_init,
    tim_file,
    flags = "",
    keys = BasicTempoKeys(),
    tparams = GeneralTempoParameter[]
    ) = BasicTempoSettings(
        work_dir,
        version,
        par_file_init,
        tim_file,
        flags,
        keys,
        tparams
        )

