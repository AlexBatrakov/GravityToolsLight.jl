#--------------------------------------------------------------------------------------------------------------
# BasicTempoSettings.jl

"""BasicTempoKeys"""
Represents the basic configuration keys for running the Tempo program.

Fields:
- `silent::Bool`: If true, Tempo will run in silent mode without printing to the standard output.
- `print_output::Bool`: If true, the output of Tempo will be printed.
- `iterative_mode::Bool`: If true, Tempo will run in an iterative mode for fitting.
- `fit_EFACs_EQUADs::Bool`: If true, Tempo will fit for EFACs and EQUAD parameters during the run.

Constructor:
- `BasicTempoKeys(;silent, print_output, iterative_mode, fit_EFACs_EQUADs)`: Creates a new `BasicTempoKeys` instance with the provided settings.
"""

```julia
struct BasicTempoKeys
    silent::Bool
    print_output::Bool
    iterative_mode::Bool
    fit_EFACs_EQUADs::Bool
    BasicTempoKeys(;silent = true, print_output = false, iterative_mode=true, fit_EFACs_EQUADs::Bool = false) = new(silent, print_output, iterative_mode, fit_EFACs_EQUADs)
end
```

"""Base.show(io::IO, keys::BasicTempoKeys)"""
Display the configuration keys for Tempo in a formatted way.

- `io::IO`: The IO stream to print to.
- `keys::BasicTempoKeys`: The keys to display.
"""

```julia
function Base.show(io::IO, keys::BasicTempoKeys)
    println(io, "Tempo keys:")
    print(io, "        Silent mode: ", keys.fit_EFACs_EQUADs)
    print(io, "        Fit EFAC and EQUAD parameters: ", keys.fit_EFACs_EQUADs)
    return nothing
end
```

#--------------------------------------------------------------------------------------------------------------
"""BasicTempoSettings{T <: AbstractTempoVersion}"""
Represents the basic settings required to run a Tempo analysis.

Fields:
- `work_dir::String`: The working directory for Tempo runs.
- `version::T`: The Tempo version to use, as a subtype of `AbstractTempoVersion`.
- `par_file_init::String`: The initial parameter file for the pulsar.
- `tim_file::String`: The file containing the times of arrival (TOAs).
- `flags::String`: Additional flags for the Tempo run.
- `keys::BasicTempoKeys`: The keys indicating the basic run settings.
- `tparams::Vector{GeneralTempoParameter}`: A vector of general Tempo parameters to be used in the run.

Constructor:
- `BasicTempoSettings(;work_dir, version, par_file_init, tim_file, flags, keys, tparams)`: Creates a new `BasicTempoSettings` instance with the provided settings.
"""

```julia
struct BasicTempoSettings{T <: AbstractTempoVersion}
    work_dir::String
    version::T
    par_file_init::String
    tim_file::String
    flags::String
    keys::BasicTempoKeys
    tparams::Vector{GeneralTempoParameter}
end
```

"""Base.show(io::IO, bsets::BasicTempoSettings)"""
Display the basic Tempo settings in a formatted way.

- `io::IO`: The IO stream to print to.
- `bsets::BasicTempoSettings`: The settings to display.
"""

```julia
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
```

```julia
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
```
