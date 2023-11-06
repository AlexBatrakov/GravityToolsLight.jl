#--------------------------------------------------------------------------------------------------------------

struct BasicTempoKeys
    silent::Bool
    print_output::Bool
    iterative_mode::Bool
    fit_EFACs_EQUADs::Bool
    BasicTempoKeys(;silent = true, print_output = false, iterative_mode=true, fit_EFACs_EQUADs::Bool = false) = new(silent, print_output, iterative_mode, fit_EFACs_EQUADs)
end

function Base.show(io::IO, keys::BasicTempoKeys)
    println(io, "Tempo keys:")
    print(io, "        Silent mode: ", keys.fit_EFACs_EQUADs)
    print(io, "        Fit EFAC and EQUAD parameters: ", keys.fit_EFACs_EQUADs)
	return nothing
end

#--------------------------------------------------------------------------------------------------------------
# Основные настройки Tempo
struct BasicTempoSettings{T <: AbstractTempoVersion}
    work_dir::String
    version::T
    par_file_init::String
    tim_file::String
    flags::String
    keys::BasicTempoKeys
    tparams::Vector{GeneralTempoParameter}
end

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

