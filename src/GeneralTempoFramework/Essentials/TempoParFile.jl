#--------------------------------------------------------------------------------------------------------------
# tempo par files

mutable struct TempoParFile
    name::String
#    path::String
    tparams::Dict{Symbol,GeneralTempoParameter}
    order::Vector{Symbol}
end

function TempoParFile()
    return TempoParFile("", Dict{Symbol,GeneralTempoParameter}(), Vector{Symbol}())
end

function TempoParFile(par_file_path::String; new_name_suffix="")
    par_file = TempoParFile(par_file_path, Dict{Symbol,GeneralTempoParameter}(), Vector{Symbol}())
    read_par_file!(par_file)
    par_file.name = par_file.name[1:end-4] * new_name_suffix * ".par"
    return par_file
end

function Base.show(io::IO, par_file::TempoParFile)
    println(io, "Tempo parameter file $(par_file.name): ")
    for (i, name_symbol) in enumerate(par_file.order)
        #print(IOContext(io, :indent => indent+4), par_file.tparams[name_symbol])
        print("    ", get_par_file_representation(par_file.tparams[name_symbol]))
        if i < length(par_file.order)
            print("\n")
        end
    end
	return nothing
end

function read_par_file!(par_file::TempoParFile)
    par_file.order = Vector{Symbol}()
    open(par_file.name, "r") do file_in
        for line in eachline(file_in)
            if startswith(line, "C ") || startswith(line, "c ")
                continue
            end
            tparam = extract_GeneralTempoParameter(line)
            par_file.tparams[tparam.name_symbol] = tparam
            push!(par_file.order, tparam.name_symbol)
        end
    end
    return par_file
end


function write_par_file(par_file::TempoParFile, name_out=par_file.name)
    open(name_out, "w") do file_out
        for name_symbol in par_file.order
            tparam = par_file.tparams[name_symbol]
            println(file_out, get_par_file_representation(tparam))
        end
    end
    return par_file
end

function modify_par_file!(par_file::TempoParFile, tparam_name::Symbol, value=tparam.value, flag=tparam.flag, uncertainty=tparam.uncertainty)
    return par_file
end

function update_par_file()

end

function extend_par_file!(par_file::TempoParFile, tparam::GeneralTempoParameter)
    if tparam.value == nothing
        par_file.tparams[tparam.name_symbol].flag = tparam.flag
    else
        if !haskey(par_file.tparams, tparam.name_symbol)
            push!(par_file.order, tparam.name_symbol)
        end
        par_file.tparams[tparam.name_symbol] = tparam
    end
    return par_file
end

function generate_par_file_path(base_name::String, suffix::String, work_dir::String)
    # Убедимся, что имя базы и суффикс не начинаются с точки и не содержат расширения файла
    base_name = strip_extension(base_name)
    suffix = lstrip(suffix, '.')

    # Сгенерируем полный путь к файлу
    if suffix != ""
        return joinpath(work_dir, base_name * "_" * suffix * ".par")
    else
        return joinpath(work_dir, base_name * ".par")
    end

    return par_file_path
end

# Убираем расширение файла, если оно есть
function strip_extension(file_name::String)
    return splitext(file_name)[1]
end

function detect_backends(tim_file_path::String)
    backends = Set{String}()  # Используем множество для избежания дубликатов
    open(tim_file_path, "r") do file
        for line in eachline(file)
            if occursin("-be ", line)
                backend_name = match(r"-be (\S+)", line)
                if backend_name !== nothing
                    # Извлекаем только название бэкенда, исключая "-be"
                    push!(backends, backend_name.captures[1])
                end
            end
        end
    end
    return collect(backends)  # Преобразуем множество в вектор
end

