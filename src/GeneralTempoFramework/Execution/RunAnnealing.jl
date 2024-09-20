function compute_energy(annealing_params::Dict{String, Float64}, basic_settings::BasicTempoSettings)
    # Формируем список параметров для Tempo2
    annealing_tparams = [TP(param_name, value, flag=0) for (param_name, value) in annealing_params]
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
    
    chisqr          = results_basic.last_internal_iteration.result.basic.fit_chisq
    rms_post_fit    = results_basic.last_internal_iteration.result.basic.rms_post_fit_residual_us
    pre_post        = results_basic.last_internal_iteration.result.basic.pre_post
    rms_tn_post_fit = results_basic.last_internal_iteration.result.basic.rms_tn_post_fit_residual_us

    println("Run finished: $(annealing_tparams):")
    println("   chisqr = $chisqr, rms_post_fit = $rms_post_fit, rms_tn_post_fit = $rms_tn_post_fit, pre_post = $pre_post")

    # Извлекаем значение хи-квадрат
    chiqsr = results_basic.last_internal_iteration.result.basic.fit_chisq
    return chiqsr
end

# Функция для вычисления энергии в параллельных процессах
function compute_energy_parallel(new_params, basic_settings, work_dir, worker_id)

    basic_settings_copy = deepcopy(basic_settings)

    # Меняем текущий каталог на каталог рабочего
    basic_settings_copy.work_dir = "$work_dir/worker$worker_id"

    # Вычисляем энергию с измененными параметрами
    return compute_energy(new_params, basic_settings_copy)
end

# Функция для инициализации графика с автоматическим извлечением параметров и их границ
function initialize_plot(annealing_settings::AnnealingSettings)
    # Извлекаем параметры для графика
    param1 = annealing_settings.parameters[1]
    param2 = annealing_settings.parameters[2]
    
    # Создаем график
    fig, ax = subplots()
    
    # Установка пределов осей по границам параметров
    ax.set_xlim(param1.min_value, param1.max_value)  # Диапазон для первого параметра (например, P_DELTA)
    ax.set_ylim(param2.min_value, param2.max_value)  # Диапазон для второго параметра (например, P_PHI)
    
    # Настройка подписей осей
    ax.set_xlabel(param1.name)  # Подпись для оси X
    ax.set_ylabel(param2.name)  # Подпись для оси Y
    
    # Пустой scatterplot
    scatterplot = ax.scatter([], [], c=[], cmap="viridis")
    
    # Добавляем цветовую шкалу
    cbar = fig.colorbar(scatterplot)
    cbar.set_label("χ² - min(χ²)")
    
    # Устанавливаем границы для цветовой шкалы
    scatterplot.set_clim(0, 10)  # Установка границ для цветовой шкалы
    
    return fig, ax, scatterplot
end

# Функция для обновления графика
function update_plot(scatterplot, x_data, y_data, energy_data, best_energy)
    # Пересчитываем цвета на основе новой энергии
    color_data = [energy - best_energy for energy in energy_data]
    
    # Обновляем позиции и цвета точек
    scatterplot.set_offsets([x_data y_data])
    scatterplot.set_array(color_data)
    
    # Перерисовываем график
    display(scatterplot.figure)
    scatterplot.figure.canvas.draw()
    scatterplot.figure.canvas.flush_events()
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

# Обновленная функция отжига с параллельными вычислениями
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
    T0 = annealing_settings.initial_temperature
    T = T0
    T_min = annealing_settings.minimum_temperature
    alpha = annealing_settings.cooling_rate
    max_iter = annealing_settings.max_iterations
    iter_per_temp = annealing_settings.iterations_per_temp

    # Инициализация текущего состояния
    num_procs = nprocs() - 1  # Количество рабочих процессов
    current_params_list = [deepcopy(Dict(param.name => param.initial_value for param in params)) for _ in 1:num_procs]
    current_energies = fill(Inf, num_procs)  # Инициализация энергий для каждого процесса
    best_params = deepcopy(current_params_list)
    best_energies = deepcopy(current_energies)

    # Автоматически определяем параметры для осей
    param1_name = params[1].name
    param2_name = params[2].name
    
    # Инициализация графика
    fig, ax, scatterplot = initialize_plot(annealing_settings)
    x_data = []  # Для хранения значений первого параметра
    y_data = []  # Для хранения значений второго параметра
    energy_data = []  # Для хранения значений энергии

    iteration = 0

    while T > T_min && iteration < max_iter
        for i in 1:iter_per_temp
            println("Iteration = $(iteration), T  = $T, T/T0 = $(T/T0), E_min = $(minimum(best_energies))")

            # Генерация новых параметров для каждого рабочего процесса
            new_params_list = Vector{Dict{String, Float64}}(undef, num_procs)
            for p in 1:num_procs
                new_params_list[p] = deepcopy(current_params_list[p])
                for param in params
                    delta_sigma = param.initial_step_size * (T / T0)
                    delta = randn() * delta_sigma
                    new_value = new_params_list[p][param.name] + delta
                    if param.is_angle
                        max_angle = param.max_value
                        new_value = mod(new_value, max_angle)
                    else
                        new_value = clamp(new_value, param.min_value, param.max_value)
                    end
                    new_params_list[p][param.name] = new_value
                end
            end

            # Использование pmap для параллельного вычисления энергии
            new_energies = pmap(p -> compute_energy_parallel(new_params_list[p], basic_settings, work_dir, p+1), 1:num_procs)
            
            # Принятие нового состояния для каждого рабочего процесса
            for (p, new_energy, new_params) in zip(1:num_procs, new_energies, new_params_list)
                delta_energy = new_energy - current_energies[p]
                if delta_energy <= 0 || rand() < exp(-delta_energy / T)
                    println("Proccess $p: E_min = $(min(best_energies[p], new_energy)). Transition! From E = $(current_energies[p]) to E = $new_energy with dE = $delta_energy.")
                    current_params_list[p] = deepcopy(new_params)
                    current_energies[p] = new_energy
                    if current_energies[p] < best_energies[p]
                        best_params[p] = deepcopy(current_params_list[p])
                        best_energies[p] = current_energies[p]
                    end
                else
                    println("Proccess $p: E_min = $(best_energies[p]). No transition!   E = $(current_energies[p]).")
                end

                # Обновление данных для графика
                push!(x_data, new_params[param1_name])  # Первый параметр
                push!(y_data, new_params[param2_name])  # Второй параметр
                push!(energy_data, new_energy)  # Сохраняем энергию
            end

            iteration += 1

            # Обновление графика с пересчетом цветов
            update_plot(scatterplot, x_data, y_data, energy_data, minimum(best_energies))

            if iteration >= max_iter
                break
            end
        end
        # Обновление температуры
        T *= alpha
    end
    return best_params, best_energies
end