#--------------------------------------------------------------------------------------------------------------
#using PyPlot

function plot_fit_results_all(output_data)
    # Получаем ключи первой итерации для определения, какие параметры имеются
    params_keys = keys(output_data[1])
    
    # Определяем количество графиков
    num_plots = length(params_keys)
    
    # Создаем сетку графиков (пример для 4 графиков в ряду)
    num_cols = 4
    num_rows = cld(num_plots, num_cols)  # округление вверх для числа строк
    
    fig, axs = subplots(num_rows, num_cols, figsize=(15, num_rows*3))
    
    # Сделаем индексацию одномерной для удобства, если всего один график
    axs = reshape(axs, num_rows*num_cols)
    
    # Итерируем по параметрам и создаем графики
    for (index, param_key) in enumerate(params_keys)
        ax = axs[index]
        values = [iteration_data[param_key] for iteration_data in output_data]
        
        # Строим график
        ax.plot(values)
        ax.set_title(param_key)
        ax.set_xlabel("Iteration")
        ax.set_ylabel("Value")
        ax.grid(true)
        
        # Если достигли последнего ключа, прекращаем отрисовку
        if index == num_plots
            break
        end
    end
    
    # Если количество графиков не заполняет всю сетку, убираем пустые subplots
    for i in (num_plots+1):(num_rows*num_cols)
        fig.delaxes(axs[i])
    end
    
    # Показываем график
    tight_layout()
    show()
end

# Теперь вызываем функцию с данными
#plot_iterations_data(output_data)

#using PyPlot

function plot_fit_results(iterations_data)
    # Извлекаем данные для графиков
    iterations = 1:length(iterations_data)
    rms_prefit = [data["RMS pre-fit residual (us)"] for data in iterations_data]
    rms_postfit = [data["RMS post-fit residual (us)"] for data in iterations_data]
    chisqr = [data["Fit Chisq"] for data in iterations_data]
    pre_post_ratio = [data["pre/post"] for data in iterations_data]
    pre_post_deviation = [ratio - 1.0 for ratio in pre_post_ratio]  # Вычитаем единицу
    
    # Создаем графики
    fig, axs = subplots(3, 1, figsize=(10, 15), constrained_layout=true)

    # # Настройка тиков оси абсцисс
    # xticks = iterations  # Убедитесь, что iterations содержит список итераций как целые числа

    # # Перебираем все подграфики и устанавливаем тики для оси X
    # for ax in axs
    #     ax.set_xticks(xticks)  # Устанавливаем тики, соответствующие итерациям
    #     ax.xaxis.set_major_formatter(plt.FuncFormatter((val, pos) -> string(Int(val)))) # Формат без дробной части
    # end

    # # Устанавливаем формат оси X для отображения целых чисел
    # for ax in axs
    #     ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    # end
    
    # График RMS pre-fit и post-fit residuals
    axs[1].semilogy(iterations, rms_prefit, label="RMS Pre-Fit")
    axs[1].semilogy(iterations, rms_postfit, label="RMS Post-Fit")
    axs[1].set_title("RMS Residuals (Log Scale)", fontsize=10)
    axs[1].legend()
    
    # График Chisqr
    axs[2].semilogy(iterations, chisqr, label="Chisqr")
    axs[2].set_title("Chisqr (Log Scale)", fontsize=10)
    axs[2].legend()
    
    # График Pre/Post Fit Ratio Deviation
    axs[3].plot(iterations, pre_post_deviation, label="Pre/Post Fit Ratio Deviation")
    axs[3].set_title("Pre/Post Fit Ratio Deviation", fontsize=10)
    axs[3].set_yscale("symlog", linthresh=1e-6)
    axs[3].legend()
    
    # Устанавливаем общие метки для осей
    for ax in axs
        ax.set_xlabel("Iteration")
        ax.set_ylabel("Value")
    end
    
    # Показываем или сохраняем фигуру
    show()
    # savefig("fit_results.png")
end

# Предполагается, что iterations_data — это массив словарей с результатами итераций
#plot_fit_results(parsed_results) # Замените parsed_results на ваш массив данных