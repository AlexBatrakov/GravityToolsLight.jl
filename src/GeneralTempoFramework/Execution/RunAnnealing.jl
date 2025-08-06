function compute_energy(annealing_params::Dict{String, Float64}, basic_settings::BasicTempoSettings)
    # Формируем список параметров для Tempo2
    annealing_tparams = [TP(param_name, value, flag=1) for (param_name, value) in annealing_params]
    # Обновляем настройки Tempo2
    annealing_basic_settings = deepcopy(basic_settings)

    # Добавляем или обновляем локальные параметры для текущей итерации
    for a_tparam in annealing_tparams
        found = false
        for (i, g_tparam) in enumerate(annealing_basic_settings.tparams)
            if g_tparam.name == a_tparam.name
                annealing_basic_settings.tparams[i] = a_tparam  # Обновляем существующий параметр
                found = true
                break
            end
        end
        if !found
            push!(annealing_basic_settings.tparams, a_tparam)  # Добавляем новый параметр
        end
    end

    println("Run started:  $(annealing_tparams):")

    # Запускаем Tempo2
    results_basic = run_tempo_basic(annealing_basic_settings)
    
    chisqr          = results_basic.last_internal_iteration.result.basic.chisqr
    rms_post_fit    = results_basic.last_internal_iteration.result.basic.rms_post_fit_residual_us
    pre_post        = results_basic.last_internal_iteration.result.basic.pre_post
    rms_tn_post_fit = results_basic.last_internal_iteration.result.basic.rms_tn_post_fit_residual_us

    println("Run finished: $(annealing_tparams):")
    println("   chisqr = $chisqr, rms_post_fit = $rms_post_fit, rms_tn_post_fit = $rms_tn_post_fit, pre_post = $pre_post")

    # Извлекаем значение хи-квадрат
    chiqsr = results_basic.last_internal_iteration.result.basic.chisqr * (1 + abs(1 - pre_post))
    return chiqsr
end

function compute_energy(annealing_params::Dict{String, Float64}, basic_settings::BasicTempoSettings, global_iters_settings::GlobalIterationsSettings)
    # Формируем список параметров для Tempo2
    annealing_tparams = [TP(param_name, value, flag=1) for (param_name, value) in annealing_params]
    # Обновляем настройки Tempo2
    annealing_basic_settings = deepcopy(basic_settings)

    # Добавляем или обновляем локальные параметры для текущей итерации
    for a_tparam in annealing_tparams
        found = false
        for (i, g_tparam) in enumerate(annealing_basic_settings.tparams)
            if g_tparam.name == a_tparam.name
                annealing_basic_settings.tparams[i] = a_tparam  # Обновляем существующий параметр
                found = true
                break
            end
        end
        if !found
            push!(annealing_basic_settings.tparams, a_tparam)  # Добавляем новый параметр
        end
    end

    println("Run started.  Init: $(annealing_params):")

    # Запускаем Tempo2
    results = run_tempo_global_iters(annealing_basic_settings, global_iters_settings)

    chisqr_init          = results.all_global_iterations[1].last_internal_iteration.result.basic.chisqr
    rms_post_fit_init    = results.all_global_iterations[1].last_internal_iteration.result.basic.rms_post_fit_residual_us
    rms_tn_post_fit_init = results.all_global_iterations[1].last_internal_iteration.result.basic.rms_tn_post_fit_residual_us
    pre_post_init        = results.all_global_iterations[1].last_internal_iteration.result.basic.pre_post

    println("Run finished. Init: $(annealing_params):")
    println("   chisqr = $chisqr_init, rms_post_fit = $rms_post_fit_init, rms_tn_post_fit = $rms_tn_post_fit_init, pre_post = $pre_post_init")
    
    chisqr_gradient           = results.last_internal_iteration.result.basic.chisqr
    rms_post_fit_gradient     = results.last_internal_iteration.result.basic.rms_post_fit_residual_us
    pre_post_gradient         = results.last_internal_iteration.result.basic.pre_post
    rms_tn_post_fit_gradient  = results.last_internal_iteration.result.basic.rms_tn_post_fit_residual_us

    annealing_params_gradient = deepcopy(annealing_params)

    try
        for param_name in keys(annealing_params_gradient)
            annealing_params_gradient[param_name] = Float64(results.final_par_file.tparams[Symbol(param_name)].value)
        end
    catch error
        println("Parfile is currupted")
    end

    println("          Gradient: $(annealing_params_gradient):")
    println("   chisqr = $chisqr_gradient, rms_post_fit = $rms_post_fit_gradient, rms_tn_post_fit = $rms_tn_post_fit_gradient, pre_post = $pre_post_gradient")
    

    E_init     = chisqr_init * (pre_post_init >= 1 ? pre_post_init : 1/pre_post_init)
    E_gradient = chisqr_gradient * (pre_post_gradient >= 1 ? pre_post_gradient : 1/pre_post_gradient)

    return E_init, E_gradient, annealing_params_gradient
end

# Функция для вычисления энергии в параллельных процессах
function compute_energy_parallel(new_params, basic_settings, work_dir, worker_id)

    basic_settings_copy = deepcopy(basic_settings)

    # Меняем текущий каталог на каталог рабочего
    basic_settings_copy.work_dir = "$work_dir/worker$worker_id"

    # Вычисляем энергию с измененными параметрами
    return compute_energy(new_params, basic_settings_copy)
end

# Функция для вычисления энергии в параллельных процессах
function compute_energy_parallel(new_params, basic_settings, global_iters_settings, work_dir, worker_id)

    basic_settings_copy = deepcopy(basic_settings)

    global_iters_settings_copy = deepcopy(global_iters_settings)

    # Меняем текущий каталог на каталог рабочего
    basic_settings_copy.work_dir = "$work_dir/worker$worker_id"

    # Вычисляем энергию с измененными параметрами
    return compute_energy(new_params, basic_settings_copy, global_iters_settings_copy)
end


# # Функция для инициализации графика с автоматическим извлечением параметров и их границ
# function initialize_plot(annealing_settings::AnnealingSettings)
#     # Извлекаем параметры для графика
#     param1 = annealing_settings.parameters[1]
#     param2 = annealing_settings.parameters[2]
    
#     # Создаем график
#     fig, ax = subplots()
    
#     # Установка пределов осей по границам параметров
#     ax.set_xlim(param1.min_value, param1.max_value)  # Диапазон для первого параметра (например, P_DELTA)
#     ax.set_ylim(param2.min_value, param2.max_value)  # Диапазон для второго параметра (например, P_PHI)
    
#     # Настройка подписей осей
#     ax.set_xlabel(param1.name)  # Подпись для оси X
#     ax.set_ylabel(param2.name)  # Подпись для оси Y
    
#     # Пустой scatterplot
#     scatterplot = ax.scatter([], [], c=[], cmap="viridis")
    
#     # Добавляем цветовую шкалу
#     cbar = fig.colorbar(scatterplot)
#     cbar.set_label("χ² - min(χ²)")
    
#     # Устанавливаем границы для цветовой шкалы
#     scatterplot.set_clim(0, 10)  # Установка границ для цветовой шкалы
    
#     return fig, ax, scatterplot
# end

# # Функция для обновления графика
# function update_plot(scatterplot, x_data, y_data, energy_data, best_energy)
#     # Пересчитываем цвета на основе новой энергии
#     color_data = [energy - best_energy for energy in energy_data]
    
#     # Обновляем позиции и цвета точек
#     scatterplot.set_offsets([x_data y_data])
#     scatterplot.set_array(color_data)
    
#     # Перерисовываем график
#     display(scatterplot.figure)
#     scatterplot.figure.canvas.draw()
#     scatterplot.figure.canvas.flush_events()
# end

# Функция для инициализации графиков для всех пар переменных
function initialize_plots(annealing_settings::AnnealingSettings)
    num_params = length(annealing_settings.parameters)  # Количество параметров
    num_plots = (num_params - 1) * (num_params - 1)  # Число графиков = (N-1) * (N-1)

    # Создаем сетку (N-1)*(N-1) под графики
    fig, axes = subplots(num_params - 1, num_params - 1, figsize=(15, 15))
    scatterplots = Dict{Tuple{Int, Int}, Any}()  # Для хранения всех scatterplot
    
    for i in 2:num_params
        for j in 1:i-1
            param1 = annealing_settings.parameters[i]
            param2 = annealing_settings.parameters[j]
            
            # Подготовка подграфика
            ax = axes[i-1, j]  # (i-1) для смещения из-за диагональной матрицы
            ax.set_xlim(param2.min_value, param2.max_value)
            ax.set_ylim(param1.min_value, param1.max_value)
            
            # Подписи только для внешних осей
            if j == 1  # Подпись слева
                ax.set_ylabel(param1.name)
            end
            if i == num_params  # Подпись снизу
                ax.set_xlabel(param2.name)
            end
            
            # Пустой scatterplot для каждой пары параметров
            scatterplot = ax.scatter([], [], c=[], cmap="viridis", alpha=0.6)
            scatterplots[(i, j)] = scatterplot  # Сохраняем scatterplot для обновлений
        end
    end

    # Удаляем ненужные оси (верхний правый угол)
    for i in 1:num_params-1
        for j in i+1:num_params-1
            axes[i, j].axis("off")  # Убираем ненужные подграфики
        end
    end
    
    # Добавляем один общий цветовой бар справа от всей сетки графиков
    fig.subplots_adjust(right=0.9)  # Оставляем место для цветовой шкалы
    cbar_ax = fig.add_axes([0.92, 0.1, 0.02, 0.8])  # Расположение для общего colorbar
    cbar = fig.colorbar(scatterplots[(2, 1)], cax=cbar_ax)  # Один colorbar для всех
    cbar.set_label("χ² - min(χ²)")
    
    # Устанавливаем начальные границы цветовой шкалы для всех графиков
    for scatterplot in values(scatterplots)
        scatterplot.set_clim(0, 100)  # Устанавливаем одинаковые границы цветовой шкалы
    end
    
    return fig, axes, scatterplots
end

# Функция для обновления всех графиков с сортировкой точек по цвету (по убыванию)
function update_plots(scatterplots, annealing_settings::AnnealingSettings, data::Dict{String,Vector{Float64}}, energy_data, best_energy)
    num_params = length(annealing_settings.parameters)
    
    # Пересчитываем цвета на основе новой энергии
    color_data = [energy - best_energy for energy in energy_data]
    
    # Сортировка данных по цвету по убыванию, чтобы точки с большими значениями цвета были сверху
    sorted_indices = sortperm(color_data, rev=true)  # Сортируем по убыванию
    sorted_color_data = color_data[sorted_indices]  # Сортируем цвета
    
    # Обновляем каждый график для всех пар параметров
    for i in 2:num_params
        for j in 1:i-1
            param1_name = annealing_settings.parameters[i].name
            param2_name = annealing_settings.parameters[j].name
            
            x_data = data[param2_name][sorted_indices]  # Сортируем x по индексу цвета
            y_data = data[param1_name][sorted_indices]  # Сортируем y по индексу цвета
            
            # Обновляем scatterplot для пары (i, j)
            scatterplot = scatterplots[(i, j)]
            scatterplot.set_offsets(hcat(x_data, y_data))  # Передаем данные x и y как матрицу
            scatterplot.set_array(sorted_color_data)  # Обновляем цвета точек
        end
    end
    
    # Перерисовываем все графики
    fig = scatterplots[(2, 1)].figure  # Получаем ссылку на фигуру
    display(fig)
    fig.canvas.draw()
    fig.canvas.flush_events()
end

# Обновленная функция отжига
function run_tempo_annealing(basic_settings::BasicTempoSettings, annealing_settings::AnnealingSettings)
    # Установка случайного зерна
    Random.seed!(annealing_settings.random_seed)

    # Извлечение настроек
    params = annealing_settings.parameters
    T = annealing_settings.initial_temperature
    T_min = annealing_settings.minimum_temperature
    alpha = annealing_settings.cooling_rate
    max_iter = annealing_settings.max_iterations
    iter_per_temp = annealing_settings.iterations_per_temp

    # Инициализация текущего состояния
    current_params = Dict(param.name => param.initial_value for param in params)
    current_energy = compute_energy(current_params, basic_settings)
    best_params = deepcopy(current_params)
    best_energy = current_energy

    # Автоматически определяем параметры для осей
    param1_name = params[1].name
    param2_name = params[2].name
    
    # Инициализация графика
    fig, ax, scatterplot = initialize_plot(annealing_settings)
    x_data = [current_params[param1_name]]  # Для хранения значений первого параметра
    y_data = [current_params[param2_name]]  # Для хранения значений второго параметра
    energy_data = [current_energy]  # Для хранения значений энергии
    update_plot(scatterplot, x_data, y_data, energy_data, best_energy)

    iteration = 0

    while T > T_min && iteration < max_iter
        for i in 1:iter_per_temp
            println("Iteration = $(iteration), temperature  = $T, best_energy = $(best_energy), current_energy = $(current_energy)")
            new_params = deepcopy(current_params)
            
            # Генерация новых параметров
            for param in params
                delta_sigma = param.initial_step_size * (T / annealing_settings.initial_temperature)
                delta = randn() * delta_sigma
                new_value = new_params[param.name] + delta
                # Обработка ограничений
                if param.is_angle
                    max_angle = param.max_value
                    new_value = mod(new_value, max_angle)
                else
                    new_value = clamp(new_value, param.min_value, param.max_value)
                end
                new_params[param.name] = new_value
            end
            
            # Вычисление новой энергии
            new_energy = compute_energy(new_params, basic_settings)
            delta_energy = new_energy - current_energy
            
            # Решение о принятии нового состояния
            if delta_energy <= 0 || rand() < exp(-delta_energy / T)
                println("Transition! From E = $current_energy to E = $new_energy with dE = $delta_energy")
                current_params = new_params
                current_energy = new_energy
                if current_energy < best_energy
                    best_params = deepcopy(current_params)
                    best_energy = current_energy
                end
            end
            iteration += 1

            # Обновление данных для графика
            push!(x_data, new_params[param1_name])  # Первый параметр
            push!(y_data, new_params[param2_name])  # Второй параметр
            push!(energy_data, new_energy)  # Сохраняем энергию
            
            # Обновление графика с пересчетом цветов
            update_plot(scatterplot, x_data, y_data, energy_data, best_energy)

            if iteration >= max_iter
                break
            end
        end
        # Обновление температуры
        T *= alpha
    end
    return best_params, best_energy
end

function run_tempo_annealing_parallel(basic_settings::BasicTempoSettings, annealing_settings::AnnealingSettings)
    work_dir = basic_settings.work_dir
    cd(work_dir)

    if annealing_settings.parallel
        for p in 2:nprocs()  # Начинаем с 2, чтобы исключить главный процесс
            rm("./worker$p", force=true, recursive=true)
            mkdir("./worker$p")
            cp(basic_settings.par_file_init, "$(work_dir)/worker$p/$(basic_settings.par_file_init)", force=true)
            cp(basic_settings.tim_file, "$(work_dir)/worker$p/$(basic_settings.tim_file)", force=true)
        end
    end

    # Установка случайного зерна
    Random.seed!(annealing_settings.random_seed)

    # Извлечение настроек
    params = annealing_settings.parameters
    energy_scale = annealing_settings.energy_scale
    Q = annealing_settings.quenching_factor
    T0 = annealing_settings.initial_temperature
    T = T0
    T_min = annealing_settings.minimum_temperature
    alpha = annealing_settings.cooling_rate
    max_iter = annealing_settings.max_iterations
    step_law = annealing_settings.step_law
    iter_per_temp = annealing_settings.iterations_per_temp

    D = length(params)
    m = log(T0 / T_min)
    n = log(max_iter)
    c = m * exp(- n * Q / D)

    # Инициализация текущего состояния
    num_procs = nprocs() - 1  # Количество рабочих процессов
    # Инициализация начальных параметров для всех процессов
    current_params_list = [
        Dict(param.name => param.min_value + rand() * (param.max_value - param.min_value) for param in params)
        for _ in 1:num_procs
    ]
    current_energies = fill(Inf, num_procs)  # Инициализация энергий для каждого процесса
    best_params = deepcopy(current_params_list)
    best_energies = deepcopy(current_energies)

    # Инициализация графиков для всех пар параметров
    fig, axes, scatterplots = initialize_plots(annealing_settings)

    # Создаем словарь для хранения значений всех параметров
    data = Dict(param.name => Float64[] for param in params)
    energy_data = []  # Для хранения значений энергии

    iteration = 0

    while T > T_min && iteration < max_iter
        for i in 1:iter_per_temp
            println("Iteration = $(iteration), T  = $T, T/T0 = $(T/T0), E_min = $(minimum(best_energies))")
            println("  scaled(T) = $(T * energy_scale), cooling rate = $(exp(-c)), k^(Q/D) = $(iteration^(Q / D))")

            # Генерация новых параметров для каждого рабочего процесса
            new_params_list = Vector{Dict{String, Float64}}(undef, num_procs)
            # Основной цикл для обновления параметров
            for p in 1:num_procs
                new_params_list[p] = deepcopy(current_params_list[p])  # Копируем текущие параметры процесса
                for param in params
                    u = rand()
                    y = sign(u - 0.5) * T * ((1 + 1/T) ^ abs(2*u - 1) - 1)

                    # Обновляем значение параметра
                    new_value = new_params_list[p][param.name] + y * (param.max_value - param.min_value)
                    new_value = update_parameter_value(param, new_value)  # Используем функцию для обработки углов и обычных параметров

                    # Присваиваем обновлённое значение параметру
                    new_params_list[p][param.name] = new_value
                end
            end

            # Использование pmap для параллельного вычисления энергии
            new_energies = pmap(p -> compute_energy_parallel(new_params_list[p], basic_settings, work_dir, p+1), 1:num_procs)
            
            # Принятие нового состояния для каждого рабочего процесса
            for (p, new_energy, new_params) in zip(1:num_procs, new_energies, new_params_list)
                delta_energy = new_energy - current_energies[p]
                if delta_energy <= 0 || rand() < exp(-delta_energy * energy_scale / T)
                    println("Worker $(p+1): E_min = $(min(best_energies[p], new_energy)). Transition! From E = $(current_energies[p]) to E = $new_energy with dE = $delta_energy.")
                    current_params_list[p] = deepcopy(new_params)
                    current_energies[p] = new_energy
                    if current_energies[p] < best_energies[p]
                        best_params[p] = deepcopy(current_params_list[p])
                        best_energies[p] = current_energies[p]
                    end
                else
                    println("Worker $(p+1): E_min = $(best_energies[p]). No transition!   E = $(current_energies[p]). Forbidden E = $new_energy")
                end

                # Обновление данных для всех параметров
                for param in params
                    push!(data[param.name], new_params[param.name])
                end
                push!(energy_data, new_energy)
            end

            iteration += 1

            # Обновление всех графиков с пересчетом цветов
            update_plots(scatterplots, annealing_settings, data, energy_data, minimum(best_energies))

            if iteration >= max_iter
                break
            end
        end
        # Обновление температуры
        T = T0 * exp(-c * iteration^(Q / D))
    end
    return best_params, best_energies
end

function run_tempo_annealing_parallel(basic_settings::BasicTempoSettings, global_iters_settings::GlobalIterationsSettings, annealing_settings::AnnealingSettings)
    work_dir = basic_settings.work_dir
    cd(work_dir)

    if annealing_settings.parallel
        for p in 2:nprocs()  # Начинаем с 2, чтобы исключить главный процесс
            rm("./worker$p", force=true, recursive=true)
            mkdir("./worker$p")
            cp(basic_settings.par_file_init, "$(work_dir)/worker$p/$(basic_settings.par_file_init)", force=true)
            cp(basic_settings.tim_file, "$(work_dir)/worker$p/$(basic_settings.tim_file)", force=true)
        end
    end

    # Установка случайного зерна
    Random.seed!(annealing_settings.random_seed)

    # Извлечение настроек
    params = annealing_settings.parameters
    energy_scale = annealing_settings.energy_scale
    Q = annealing_settings.quenching_factor
    T0 = annealing_settings.initial_temperature
    T = T0
    T_min = annealing_settings.minimum_temperature
    alpha = annealing_settings.cooling_rate
    max_iter = annealing_settings.max_iterations
    step_law = annealing_settings.step_law
    iter_per_temp = annealing_settings.iterations_per_temp

    D = length(params)
    m = log(T0 / T_min)
    n = log(max_iter)
    c = m * exp(- n * Q / D)

    # Инициализация текущего состояния
    num_procs = nprocs() - 1  # Количество рабочих процессов
    # Инициализация начальных параметров для всех процессов
    current_params_list = [
        Dict(param.name => param.min_value + rand() * (param.max_value - param.min_value) for param in params)
        for _ in 1:num_procs
    ]
    current_energies = fill(Inf, num_procs)  # Инициализация энергий для каждого процесса
    best_params = deepcopy(current_params_list)
    best_energies = deepcopy(current_energies)

    # Инициализация графиков для всех пар параметров
    fig, axes, scatterplots = initialize_plots(annealing_settings)

    # Создаем словарь для хранения значений всех параметров
    data = Dict(param.name => Float64[] for param in params)
    energy_data = []  # Для хранения значений энергии

    iteration = 0

    while T > T_min && iteration < max_iter
        for i in 1:iter_per_temp
            println("Iteration = $(iteration), T  = $T, T/T0 = $(T/T0), E_min = $(minimum(best_energies))")
            println("  scaled(T) = $(T * energy_scale), cooling rate = $(exp(-c)), k^(Q/D) = $(iteration^(Q / D))")

            # Генерация новых параметров для каждого рабочего процесса
            new_params_list_init = Vector{Dict{String, Float64}}(undef, num_procs)
            # Основной цикл для обновления параметров
            for p in 1:num_procs
                new_params_list_init[p] = deepcopy(current_params_list[p])  # Копируем текущие параметры процесса
                for param in params
                    new_value = Inf
                    while new_value < param.min_value || new_value > param.max_value
                        u = rand()
                        y = sign(u - 0.5) * T * ((1 + 1/T) ^ abs(2*u - 1) - 1)

                        # Обновляем значение параметра
                        new_value = new_params_list_init[p][param.name] + y * (param.max_value - param.min_value)
                        # new_value = update_parameter_value(param, new_value)  # Используем функцию для обработки углов и обычных параметров
                    end
                    # Присваиваем обновлённое значение параметру
                    new_params_list_init[p][param.name] = new_value
                end
            end

            # Использование pmap для параллельного вычисления энергии
            results = pmap(p -> compute_energy_parallel(new_params_list_init[p], basic_settings, global_iters_settings, work_dir, p+1), 1:num_procs)

            new_params_list_gradient = deepcopy(new_params_list_init)
            new_energies_init        = Vector{Float64}(undef, num_procs)
            new_energies_gradient    = Vector{Float64}(undef, num_procs)
            for p in 1:num_procs
                new_energies_init[p] = results[p][1]
                new_energies_gradient[p]      = results[p][2]
                if annealing_settings.gradient
                    new_params_list_gradient[p] = results[p][3]
                    for param in params
                        new_params_list_gradient[p][param.name] = update_parameter_value(param, new_params_list_gradient[p][param.name])
                    end
                end
            end

            
            # Принятие нового состояния для каждого рабочего процесса
            for (p, new_energy_init, new_params_init, new_energy_gradient, new_params_gradient) in zip(1:num_procs, new_energies_init, new_params_list_init, new_energies_gradient, new_params_list_gradient)
                for (new_energy, new_params) in [(new_energy_init, new_params_init), (new_energy_gradient, new_params_gradient)]
                    delta_energy = new_energy - current_energies[p]
                    if delta_energy <= 0 || rand() < exp(-delta_energy / (T * energy_scale))
                        println("Worker $(p+1): E_min = $(min(best_energies[p], new_energy)). Transition! From E = $(current_energies[p]) to E = $new_energy with dE = $delta_energy.")
                        current_params_list[p] = deepcopy(new_params)
                        current_energies[p] = new_energy
                        if current_energies[p] < best_energies[p]
                            best_params[p] = deepcopy(current_params_list[p])
                            best_energies[p] = current_energies[p]
                        end
                    else
                        println("Worker $(p+1): E_min = $(best_energies[p]). No transition!   E = $(current_energies[p]). Forbidden E = $new_energy")
                    end
                end
                # Обновление данных для всех параметров
                for param in params
                    push!(data[param.name], new_params_init[param.name])
                    push!(data[param.name], new_params_gradient[param.name])
                end
                push!(energy_data, new_energy_init)
                push!(energy_data, new_energy_gradient)
            end

            iteration += 1

            # Обновление всех графиков с пересчетом цветов
            update_plots(scatterplots, annealing_settings, data, energy_data, minimum(best_energies))

            if iteration >= max_iter
                break
            end
        end
        # Обновление температуры
        T = T0 * exp(-c * iteration^(Q / D))
    end
    return best_params, best_energies, scatterplots, data, energy_data
end