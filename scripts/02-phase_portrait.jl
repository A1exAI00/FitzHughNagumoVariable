#=
В систему ФитцХью-Нагумо входят параметры a,b,ε.
Пусть параметры a и ε являются постоянными, можем изменять параметр b.

В зависимости от значения параметра b будет изменяться форма предельного цикла.

В данном скрипте генерируется gif-файл, состоящий из фазовых портретов 
системы ФитцХью-Нагумо при различных значениях параметра b.
На фазовых портретах изображены кривые медленных движений, кривые быстрых движений, 
траектории системы ФитцХью-Нагумо, полученная численным интегрированием.
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
PLOT_FILENAME = "02-phase_space_fr"
PLOT_PX_PER_UNIT_PNG = 2

#########################################################################################

# Постоянные параметры системы
a, ε = 0.5, 0.001
b_min, b_max, N_b = 0.500001, 0.7, 150

# Начальные условия системы
x₀, y₀ = 0.0, 0.0

# Время интегрирования не постоянное, разные значения ε -> разные периоды 
t_span = [0.0,1e4]

#########################################################################################

b_range = range(b_min, b_max, N_b)
U₀ = [x₀, y₀]

#########################################################################################

fig = Figure(size=PLOT_RES)
ax = beautiful_axis(fig[1,1], title="Phase space", xlabel="x", ylabel="y")

record(fig, "02-phase_space.gif", eachindex(b_range); framerate=30) do i
    println(i)

    b = b_range[i]
    sys = FitzHugh_Nagumo_base_system((a,b,ε), U₀)
    system_integrate!(sys, t_span)
    sol = sys.sol

    xA_, xB_, xC_, xD_, yB_, yD_ = xA(a,b), xB(a,b), xC(a,b), xD(a,b), yB(a,b), yD(a,b)

    empty!(ax)
    ax.title = "Phase space; ε=$(ε), b=$(b)"
    limits!(ax, xB_*1.1, xD_*1.1, yD_*1.1, yB_*1.1)

    hlines!(ax, 0.0, color=:black)
    vlines!(ax, 0.0, color=:black)

    # График аналитического медленного движения
    x_slow = range(xB_, xD_, 100)
    lines!(ax, x_slow, f_slow.(x_slow,a,b), color=:red)

    # График аналитического быстрого движения
    x_fast_upper = range(xB_, xA_, 100)
    lines!(ax, x_fast_upper, 0 .*x_fast_upper .+ yB_, color=:red)
    x_fast_lower = range(xC_, xD_, 100)
    lines!(ax, x_fast_lower, 0 .*x_fast_lower .+ yD_, color=:red)

    # График предельного цикла модели
    lines!(ax, sol[1,:], sol[2,:], color=:blue)

    savingpath = joinpath(PLOT_SAVING_DIR, PLOT_FILENAME*"$(lpad(i,4)).png")
    # save(savingpath, fig, px_per_unit=PLOT_PX_PER_UNIT_PNG)
end