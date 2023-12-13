#=
В систему ФитцХью-Нагумо входят параметры a,b,ε.
Пусть параметры a и ε являются постоянными, можем изменять параметр b.

Сделаем параметр b новой фазовой переменной. 
Дифференциальное уравнение, описывающее эволюцию b, зададим уравнением 
первого порядка на прямой с 1 устойчивой точкой b=k с функцией активации 
(то есть до некоторого времени t=t_crit эволюции b не будет, а после эта 
фазовая переменная начнет эволюционировать).

Функцию активации аппроксимировал функцией гиперболического тангенса для избежания разрыва.
=#

#########################################################################################

using DrWatson
@quickactivate "FitzHughNagumoVariable"

using CairoMakie

include(srcdir("FitzHughNagumoVariable.jl"))
using .FitzHughNagumoVariable

#########################################################################################

# Настройки генерируемого графика
PLOT_RES = (1000, 800)
PLOT_SAVING_DIR = plotsdir(); println(pwd())
PLOT_FILENAME = "04-system1-"
savingpath = joinpath(PLOT_SAVING_DIR, PLOT_FILENAME*"$(time_ns()).png")
PLOT_PX_PER_UNIT_PNG = 2

#########################################################################################

# Постоянные параметры системы
a, ε  = 0.5, 1e-3
μ, k = 1e-1, 2.0

# Начальные условия системы
x₀, y₀ = 0.0, 0.0
b₀_min, b₀_max, N_b₀ = 0.6, 3.0, 10

# Время интегрирования
t_min = 0
t_max = 100 * T_Σ_analitic(a, b₀_max, ε)
t_crit = (t_min + t_max)/2

#########################################################################################

b₀_range = range(b₀_min, b₀_max, N_b₀)
system_param = (a, ε, μ, k, t_crit)
t_span = [t_min, t_max]

#########################################################################################

T_s = Any[]
t_sols = Any[]
for i in eachindex(b₀_range)
    println("i = $(i)")

    # Выбрать b₀ для начальных условий
    b₀ = b₀_range[i]

    sys = FitzHugh_Nagumo_1_system(system_param, [x₀, y₀, b₀])
    system_integrate!(sys, t_span)
    sol = sys.sol

    x_sol = sol[1,:]
    t_sol = sol.t

    push!(T_s, measure_period(x_sol, t_sol))
    push!(t_sols, t_sol)
end

#########################################################################################

fig = Figure(resolution=PLOT_RES)
ax = Axis(fig[1,1], 
    title="T(t),",
    xlabel="t",
    ylabel="T",
    xminorticksvisible = true, 
	xminorgridvisible = true, 
	yminorticksvisible = true, 
	yminorgridvisible = true, 
    xminorticks = IntervalsBetween(10),
	yminorticks = IntervalsBetween(10)
)

hlines!(ax, 0.0, color=:black)
vlines!(ax, 0.0, color=:black)

for i in eachindex(b₀_range)
    lines!(ax, t_sols[i], T_s[i], label="b₀=$(b₀_range[i])")
end

N_b₀ ≤ 10 ? axislegend(ax, position=:rt) : nothing
# if N < 7 axislegend(ax, position=:rt) end
save(savingpath, fig, px_per_unit=PLOT_PX_PER_UNIT_PNG)