#=
В систему ФитцХью-Нагумо входят параметры a,b,ε.
Пусть параметры a и ε являются постоянными, можем изменять параметр b.

При условии малосьи ε можно аналитически рассчитать период 
предельного цикла в системе ФитцХью-Нагумо в зависимости от параметров a,b,ε 
(в приближении того, что время прохождения быстрых движений равно 0).

В данном скрипте сравниваются периоды, посчитанные аналитически, и полученные 
путем численного интегрирования системы ФитцХью-Нагумо.

Показываются 1 графика (2 файла в разных масштабах осей):
Значение периода в зависимости от значения параметра T(b) в диапозоне b∈[b_min, b_max] 
для набора значений ε∈[ε_min, ε_max].
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
# TODO refactor naming and saving names
PLOT_RES = (1000, 800)
PLOT_SAVING_DIR = plotsdir(); println(pwd())
PLOT_FILENAME_1 = "01-T(b)_lin"
PLOT_FILENAME_2 = "01-T(b)_log"
savingpath_1 = joinpath(PLOT_SAVING_DIR, PLOT_FILENAME_1*"$(time_ns()).png")
savingpath_2 = joinpath(PLOT_SAVING_DIR, PLOT_FILENAME_2*"$(time_ns()).png")
PLOT_PX_PER_UNIT_PNG = 2

#########################################################################################

# Постоянные параметры системы
a = 0.5
ε_min, ε_max, N_ε = 1e-4, 1.0, 5
b_min, b_max, N_b = a+0.1, 7.0, 30

# Начальные условия системы
x₀, y₀ = 0, 0

# Время интегрирования не постоянное, разные значения ε -> разные периоды 
t_model_mult = 50

#########################################################################################

ε_range = exp.(range(log(ε_min), log(ε_max), N_ε))
b_range = exp.(range(log(b_min), log(b_max), N_b))
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
ax = beautiful_axis(fig[1,1], title="T(b), linear scale", xlabel="b", ylabel="T")

hlines!(ax, 0.0, color=:black)

for i in eachindex(ε_range)
    lines!(ax, b_range, T_analitic[:,i], label="ε=$(round(ε_range[i], digits=5))")
    scatter!(ax, b_range, T_model[:,i])
end

axislegend(ax, position=:lt)
save(savingpath_1, fig, px_per_unit=PLOT_PX_PER_UNIT_PNG)

#########################################################################################

fig = Figure(size=PLOT_RES)
ax = beautiful_axis(fig[1,1], title="T(b), log-log scale", xlabel="b", ylabel="T", is_x_log=true, is_y_log=true)

for i in eachindex(ε_range)
    lines!(ax, b_range, T_analitic[:,i], label="ε=$(round(ε_range[i], digits=5))")
    scatter!(ax, b_range, T_model[:,i])
end

axislegend(ax, position=:lt)
save(savingpath_2, fig, px_per_unit=PLOT_PX_PER_UNIT_PNG)