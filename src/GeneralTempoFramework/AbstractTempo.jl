#--------------------------------------------------------------------------------------------------------------
# Определение абстрактного типа для версий Tempo
abstract type AbstractTempoVersion end

# Определение функции для всех типов Tempo
function get_tempo_command(::AbstractTempoVersion)
    error("This function must be implemented for each specific Tempo version.")
end

# Определение конструкторов и функций для Tempo
struct Tempo <: AbstractTempoVersion
    custom_data_directory::String
    Tempo() = new(get_tempo_directory()) # Конструктор по умолчанию
    Tempo(custom_dir::String) = new(custom_dir) # Конструктор с параметром
end

# Функция для получения пути к директории по умолчанию для Tempo
get_tempo_directory() = get(ENV, "TEMPO", "/default/path/to/tempo")

# Определение конструкторов и функций для Tempo2
struct Tempo2 <: AbstractTempoVersion
    custom_data_directory::String
    Tempo2() = new(get_tempo2_directory()) # Конструктор по умолчанию
    Tempo2(custom_dir::String) = new(custom_dir) # Конструктор с параметром
end

# Функция для получения пути к директории по умолчанию для Tempo2
get_tempo2_directory() = get(ENV, "TEMPO2", "/default/path/to/tempo2")

# Определение функций для получения команды
get_tempo_command(::Tempo) = "tempo"
get_tempo_command(::Tempo2) = "tempo2"

