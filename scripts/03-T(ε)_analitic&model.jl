#=
В систему ФитцХью-Нагумо входят параметры a,b,ε.
Пусть параметры a и b являются постоянными, можем изменять параметр ε.

При условии малосьи ε можно аналитически рассчитать период 
предельного цикла в системе ФитцХью-Нагумо в зависимости от параметров a,b,ε 
(в приближении того, что время прохождения быстрых движений равно 0).

В данном скрипте сравниваются периоды, посчитанные аналитически, и полученные 
путем численного интегрирования системы ФитцХью-Нагумо.

Показываются 1 графика (2 файла в разных масштабах осей):
Значение периода в зависимости от значения параметра T(ε) в диапозоне ε∈[ε_min, ε_max]
для набора значений b∈[b_min, b_max].
Кружочками отмечены результаты численного интегрирования. Линией отмечены аналитические результаты.
=#

#########################################################################################

using DrWatson
@quickactivate "FitzHughNagumoVariable"

include(srcdir("plots.jl"))

include(srcdir("FitzHughNagumoVariable.jl"))
using .FitzHughNagumoVariable

#########################################################################################

# Настройки генерируемого графика
PLOT_RES = (1000, 800)
PLOT_FILENAME_1 = "03-T(ε)_lin$(time_ns()).png"
PLOT_FILENAME_2 = "03-T(ε)_log$(time_ns()).png"
PLOT_SAVE_PATH_1 = plotsdir(PLOT_FILENAME_1)
PLOT_SAVE_PATH_2 = plotsdir(PLOT_FILENAME_2)
PLOT_PX_PER_UNIT_PNG = 2

#########################################################################################

# Постоянные параметры системы
a = 0.5
ε_min, ε_max, N_ε = 1e-6, 1e-3, 20
b_min, b_max, N_b = 1, 4, 10

# Начальные условия системы
x₀, y₀ = 0.0, 0.0

# Время интегрирования не постоянное, разные значения ε -> разные периоды
t_model_mult = 50

#########################################################################################

ε_range = exp.(range(log(ε_min), log(ε_max), N_ε))
b_range = range(b_min, b_max, N_b)
U₀ = [x₀, y₀]

#########################################################################################

# Аналитический период
T_analitic = [T_Σ_analitic(a,b,ε) for b in b_range, ε in ε_range]

T_model = zeros((N_b, N_ε))
for i in eachindex(b_range), j in eachindex(ε_range)
    println(i, " ", j)
    
    # Выбор параметров из наборов
    b, ε = b_range[i], ε_range[j]
    # Переменное время интегрирования, чтобы поместилось много периодов
    t_span = [0, t_model_mult*T_analitic[i,j]] 

    sys = FitzHugh_Nagumo_base_system((a,b,ε), U₀)
    system_integrate!(sys, t_span)
    sol = sys.sol

    T_model[i,j] = measure_period_avg(sol[1,:], sol.t)
end

#########################################################################################

fig = Figure(size=PLOT_RES)
ax = beautiful_axis(fig[1,1], title="T(ε), linear scale", xlabel="ε", ylabel="T")

hlines!(ax, 0.0, color=:black)
hlines!(ax, 0.0, color=:black)

for i in eachindex(b_range) 
    lines!(ax, ε_range, T_analitic[i,:], label="b=$(b_range[i])")
    scatter!(ax, ε_range, T_model[i,:])
end

axislegend(ax)
save(PLOT_SAVE_PATH_1, fig, px_per_unit=PLOT_PX_PER_UNIT_PNG)

#########################################################################################

fig = Figure(size=PLOT_RES)
ax = beautiful_axis(fig[1,1], title="T(ε), log-log scale", xlabel="ε", ylabel="T", is_x_log=true, is_y_log=true)

for i in eachindex(b_range) 
    lines!(ax, ε_range, T_analitic[i,:], label="b=$(b_range[i])")
    scatter!(ax, ε_range, T_model[i,:])
end

axislegend(ax)
save(PLOT_SAVE_PATH_2, fig, px_per_unit=PLOT_PX_PER_UNIT_PNG)