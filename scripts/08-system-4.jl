#=
В данном скрипте реализовано то же самое, что 05-model2_learning.jl, 
но разделяя систему на N+1: эволюция k во времени и N раздельных осцилляторов.
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
PLOT_RES = (1000, 1600)
PLOT_SAVING_DIR = plotsdir(); println(pwd())
PLOT_FILENAME = "08-system-4-"
savingpath = joinpath(PLOT_SAVING_DIR, PLOT_FILENAME*"$(time_ns()).png")
PLOT_PX_PER_UNIT_PNG = 2

#########################################################################################

# Постоянные параметры системы
N = 100 # Количество осцилляторов
a = 0.5
ε = 1e-3
μ = 1e-1
d = 10.0 # Коэффициент связи между элементами цепочки
Δt = 100.0 # Временной масштаб функции активации
f(x) = (x-1.0)*(2.0-x)*(x-1.5)

# Начальные условия системы
x₀, y₀, b₀ = 0.0, 0.0, 1.5
k₀_min, k₀_max = 0.6, 2.51

# Время интегрирования
t_min = 0
t_max = 50 * T_Σ_analitic(a, b₀, ε)
t₀ = (t_max-t_min)/10

#########################################################################################

aⱼ = a .* ones(N)
system_param = (N, t₀, Δt, ε, μ, d, f, aⱼ...)
x₀ⱼ = x₀ .* ones(N)
y₀ⱼ = y₀ .* ones(N)
b₀ⱼ = b₀ .* ones(N)
k₀ⱼ = range(k₀_min, k₀_max, N)
U₀ = [x₀ⱼ..., y₀ⱼ..., b₀ⱼ...,k₀ⱼ...]
t_span = (t_min, t_max)

println("system_param=$(system_param)")
println("U₀=$(U₀)")

#########################################################################################

sys = FitzHugh_Nagumo_4_system(system_param, U₀)
system_integrate!(sys, t_span)
sols = sys.sols
k_sol = sys.k_sol

T = [measure_period(sols[i][1,:], sols[i].t) for i in 1:N]

#########################################################################################

fig = Figure(size=PLOT_RES)
ax_x = beautiful_axis(fig[1,1], title="x(t)", xlabel="t", ylabel="x")
ax_b = beautiful_axis(fig[2,1], title="b(t)", xlabel="t", ylabel="b")
ax_k = beautiful_axis(fig[3,1], title="k(t), d=$(d)", xlabel="t", ylabel="k")
ax_T = beautiful_axis(fig[4,1], title="T(t)", xlabel="t", ylabel="T")
ax_b_n = Axis(fig[5,1], title="b(n), d=$(d)", xlabel="n", ylabel="b")

# axis
hlines!.((ax_x, ax_b, ax_k, ax_T, ax_b_n), 0.0, color=:black)
vlines!.((ax_x, ax_b, ax_k, ax_T, ax_b_n), 0.0, color=:black)

for i in 1:N
    lines!(ax_x, sols[i].t, sols[i][1,:], label="$(i)")
    lines!(ax_b, sols[i].t, sols[i][3,:], label="")
    lines!(ax_k, k_sol.t, k_sol[i,:], label="")
    lines!(ax_T, sols[i].t, T[i], label="")
end
scatter!(ax_b_n, 1:N, k_sol[:,end])

# TODO refactor
if N < 10
    axislegend(ax_x, position=:rt)
end

save(savingpath, fig, px_per_unit=PLOT_PX_PER_UNIT_PNG)