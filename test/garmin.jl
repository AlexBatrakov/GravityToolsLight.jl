using CSV
using DataFrames
using PyPlot
using Dates
using PyCall

@pyimport matplotlib.dates as mdates
@pyimport matplotlib.ticker as ticker

# Функция для конвертации времени из формата "Xh Ymin" в минуты
function duration_to_minutes(duration_str)
    if duration_str == "--"
        return nothing
    end
    duration_str = replace(duration_str, "min" => "", "m" => "")
    parts = split(duration_str, "h ")
    hours = 0
    minutes = 0
    if length(parts) == 2
        hours, minutes = parse.(Int, parts)
    elseif occursin("h", duration_str)
        hours = parse(Int, parts[1])
    else
        minutes = parse(Int, parts[1])
    end
    return hours * 60 + minutes
end


cd("/Users/abatrakov/Documents/FUN/garmin_data")

# Чтение и обработка данных о стрессе
stress_data = CSV.read("stress.dat", DataFrame, delim='\t', header=false)
rename!(stress_data, [:date, :stress, :rest, :low, :medium, :high])
stress_data.rest = duration_to_minutes.(stress_data.rest)
stress_data.low = duration_to_minutes.(stress_data.low)
stress_data.medium = duration_to_minutes.(stress_data.medium)
stress_data.high = duration_to_minutes.(stress_data.high)

# Чтение и обработка данных о сне
sleep_data = CSV.read("sleep.dat", DataFrame, delim='\t', header=false)
rename!(sleep_data, [:date, :score, :quality, :duration, :bedtime, :wake_time])
sleep_data.duration = duration_to_minutes.(sleep_data.duration)

# Чтение и обработка данных о пульсе покоя
resting_heart_rate_data = CSV.read("resting_heart_rate.dat", DataFrame, delim='\t', header=false)
rename!(resting_heart_rate_data, [:date, :resting_heart_rate])
resting_heart_rate_data.resting_heart_rate = parse.(Int, replace.(resting_heart_rate_data.resting_heart_rate, " bpm" => ""))

# Вывод первых строк данных для проверки
first(stress_data, 5), first(sleep_data, 5), first(resting_heart_rate_data, 5)


# Обновленная функция для проверки и преобразования дат
function check_and_convert_dates(df::DataFrame, date_col::Symbol)
    # Проверяем, является ли первый элемент столбца даты объектом Date или DateTime
    if typeof(df[1, date_col]) <: Date
        return df[!, date_col]  # Если да, ничего не делаем
    else
        # Попытка преобразования строк в объекты Date
        try
            return Date.(df[!, date_col], "yyyy-mm-dd")
        catch e
            println("Ошибка при преобразовании дат: ", e)
            return df[!, date_col]  # Возвращаем исходные данные в случае ошибки
        end
    end
end

# Применяем функцию к каждому DataFrame
stress_data[!, :date] = check_and_convert_dates(stress_data, :date)
sleep_data[!, :date] = check_and_convert_dates(sleep_data, :date)
resting_heart_rate_data[!, :date] = check_and_convert_dates(resting_heart_rate_data, :date)

# Проверяем, что даты преобразованы корректно
println(stress_data.date)
println(sleep_data.date)
println(resting_heart_rate_data.date)

# Предположим, что столбец :stress в stress_data содержит уровни стресса
# И преобразуем даты в числовой формат для PyPlot, если это еще не сделано
stress_levels = stress_data.stress # замените :stress на актуальное название столбца, если оно отличается
sleep_scores =  sleep_data.score
convert_to_float_or_nan = x -> x != "--" ? parse(Float64, x) : NaN
sleep_scores = map(convert_to_float_or_nan, sleep_scores)
resting_heart_rates = resting_heart_rate_data.resting_heart_rate









# Построение графиков
fig, axs = plt.subplots(3, 1, figsize=(10, 8)) # создаем 3 подграфика

# График стресса
axs[1].plot(stress_data.date, stress_levels, "r-")
axs[1].xaxis.set_major_locator(mdates.AutoDateLocator(interval=10)) # каждые 10 дней
axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
axs[1].grid(true)

# График качества сна
axs[2].bar(sleep_data.date, sleep_scores, color="C1")
axs[2].set_yticks(0:20:100) # Тики от 0 до 100 с шагом 20
axs[2].xaxis.set_major_locator(ticker.MaxNLocator(nbins=5)) # Ограничение количества тиков
axs[2].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
axs[2].grid(true)

# График пульса покоя
axs[3].plot(resting_heart_rate_data.date, resting_heart_rates, "b-")
axs[3].xaxis.set_major_locator(ticker.MaxNLocator(nbins=5)) # Ограничение количества тиков
axs[3].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
axs[3].grid(true)


# Проверяем, что даты преобразованы в числовой формат
stress_data.date = mdates.date2num.(stress_data.date)

# Настройка тиков на оси X для второго подграфика (качество сна)
axs[2].xaxis.set_major_locator(mdates.AutoDateLocator())
axs[2].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
axs[2].locator_params(axis='x', nbins=5) # Ограничим количество тиков

# Настройка тиков на оси X для третьего подграфика (пульс покоя)
axs[3].xaxis.set_major_locator(mdates.AutoDateLocator())
axs[3].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
axs[3].locator_params(axis='x', nbins=5) # Ограничим количество тиков

plt.tight_layout()
plt.show()
















# Объединяем все даты и убираем дубликаты
all_dates = unique(vcat(stress_data.date, sleep_data.date, resting_heart_rate_data.date))

# Сортируем даты
sort!(all_dates)

# Функция для получения значения с проверкой на наличие данных
function get_value_or_nan(date, data, column)
    row = findfirst(==(date), data.date)
    if isnothing(row)
        return NaN
    else
        value = data[row, column]
        # Check if the value is a string and convert it to a number or NaN
        if typeof(value) == String
            return tryparse(Float64, value) === nothing ? NaN : parse(Float64, value)
        else
            return value
        end
    end
end

# Создаём массивы данных с значениями или NaN
dates = String[]
stress_levels = Float64[]
sleep_scores = Float64[]
resting_heart_rates = Float64[]

for date in all_dates
    push!(dates, Dates.format(date, "Y-m-d"))
    push!(stress_levels, get_value_or_nan(date, stress_data, :stress))
    push!(sleep_scores, get_value_or_nan(date, sleep_data, :score))
    push!(resting_heart_rates, get_value_or_nan(date, resting_heart_rate_data, :resting_heart_rate))
end

# Теперь массивы dates, stress_levels, sleep_scores и resting_heart_rates имеют одинаковую длину
# и содержат NaN вместо отсутствующих значений


















# Функция для преобразования строки даты в объект Date в Julia
function parse_date(date_str)
    try
        return Date(date_str, "d u")
    catch
        return nothing # В случае ошибки парсинга вернуть nothing
    end
end

# Преобразование строк с датами в объекты Date для каждого DataFrame
stress_data.date = parse_date.(stress_data.date)
sleep_data.date = parse_date.(sleep_data.date)
resting_heart_rate_data.date = parse_date.(resting_heart_rate_data.date)

# Удаление записей, где дата является `nothing`
stress_data = stress_data[.!isnothing.(stress_data.date), :]
sleep_data = sleep_data[.!isnothing.(sleep_data.date), :]
resting_heart_rate_data = resting_heart_rate_data[.!isnothing.(resting_heart_rate_data.date), :]

# Подготовим данные, исключив `nothing`
valid_indices = .!isnothing.(stress_data.date)
dates = [Dates.format(d, "Y-m-d") for d in stress_data.date[valid_indices]]
stress_levels = stress_data.stress[valid_indices]

# Создаем новый график
figure(figsize=(10, 4))  # Указываем размер фигуры для улучшения видимости

# Построение графика
plot_date(dates, stress_levels, "-")


# Настройка меток и заголовка
xlabel("Date")
ylabel("Stress Level")
title("Stress Level Over Time")

# Поворачиваем тики на оси X
xticks(rotation=45)

# Делаем тики реже, если их слишком много
if length(dates) > 30
    # Выбираем каждый n-ый тик для отображения
    n = ceil(Int, length(dates) / 10)  # для примера, где n - это интервал между тиками
    xticks(1:n:length(dates), dates[1:n:end])
end

tight_layout()  # Автоматически подгоняет размеры элементов, чтобы они не перекрывались

# Отображаем график
show()

# ... (начало скрипта остается без изменений)
# ... (начало скрипта остается без изменений)

# Создаем окно с несколькими подграфиками
fig, axs = subplots(3, 1, figsize=(12, 10))  # 3 подграфика

# Устанавливаем общий заголовок для всех графиков
suptitle("Health Data Overview")

# График уровня стресса
axs[1].plot_date(dates, stress_levels, "-", color="C0")
axs[1].set_ylabel("Stress Level")
axs[1].tick_params(axis="y", labelcolor="C0")
axs[1].set_xticklabels([])  # Убрать метки на оси X для этого подграфика

# График качества сна
sleep_quality_levels = [sleep_data.score[valid_indices]...]  # предполагаем, что это числовые значения
axs[2].bar(dates, sleep_quality_levels, color="C1")
axs[2].set_ylabel("Sleep Quality Score")
axs[2].set_xticklabels([])  # Убрать метки на оси X для этого подграфика

# График пульса в покое
# Предполагаем, что resting_heart_rate_data уже отфильтрованы по датам
resting_heart_rate_levels = resting_heart_rate_data.resting_heart_rate
axs[3].plot_date(dates, resting_heart_rate_levels, "-", color="C2")
axs[3].set_ylabel("Resting Heart Rate (bpm)")
axs[3].tick_params(axis="y", labelcolor="C2")

# Форматирование оси X для нижнего графика
# Устанавливаем формат даты для каждого подграфика
date_formatter = mdates.DateFormatter("%Y-%m-%d")
for ax in axs
    ax.xaxis[:set_major_formatter](date_formatter)
end

# Регулируем тики, чтобы они не были слишком густыми
if length(dates) > 30
    step = ceil(Int, length(dates) / 10)
    for ax in axs
        ax.set_xticks(ax.get_xticks()[1:step:end])
    end
end

# Поворот меток на оси X
for label in axs[3].get_xticklabels()
    label.set_rotation(45)
    label.set_ha("right")  # горизонтальное выравнивание
end

# Автоматический подгон размера графиков, чтобы надписи не обрезались
fig.autofmt_xdate()
fig.tight_layout(rect=[0, 0.03, 1, 0.95])  # rect - чтобы оставить место для общего заголовка

# Отображаем график
show()






