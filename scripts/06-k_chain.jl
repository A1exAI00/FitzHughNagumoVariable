#=
Из параграфа Chaotic and Regular Patterns in Two-Dimentional Lattices of Coupled Bistable Units
была взята модель бистабильного элемента kⱼ.
В ходе создания цепочки из таких элементов они приходят к разным состояниям равновесия, 
и в зависимости от величины их связи d можно добиться сильного "расщипления". 

В данном скрипте интегрируется цепочка из N бистабильных элементов.
Показываются 2 графика:
1. Эпюры kⱼ(t)
2. Распределеление kⱼ по цепочке в последний момент времени после интегрирования
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
PLOT_FILENAME = "06-k_chain-"
savingpath = joinpath(PLOT_SAVING_DIR, PLOT_FILENAME*"$(time_ns()).png")
PLOT_PX_PER_UNIT_PNG = 2

#########################################################################################

# Постоянные параметры системы
N, ε, μ, a, d = 100, 1e-3, 1e-1, 0.5, 30.0
step_time_scale = 100
f(x) = (x-1.0)*(2.0-x)*(x-1.5)

# Начальные условия системы
x₀, y₀, b₀ = 0.0, 0.0, 1.5
k₀_min, k₀_max = 0.6, 2.5

# Время интегрирования
t_min = 0
t_max = 50 * T_Σ_analitic(a, b₀, ε)
t₀ = (t_max-t_min)/10

#########################################################################################

aⱼ = a .* ones(N)
system_param = (N, t₀, step_time_scale, ε, μ, d, f, aⱼ...)
x₀ⱼ = x₀ .* ones(N)
y₀ⱼ = y₀ .* ones(N)
b₀ⱼ = b₀ .* ones(N)
k₀ⱼ = range(k₀_min, k₀_max, N)
U₀ = [x₀ⱼ..., y₀ⱼ..., b₀ⱼ...,k₀ⱼ...]
t_span = (t_min, t_max)

@show (N, d, f)

#########################################################################################

sys = FitzHugh_Nagumo_3_system(system_param, U₀, false)
system_integrate!(sys, t_span)
k_sol = sys.k_sol
t_sol = k_sol.t

k_final = k_sol[:,end]

#########################################################################################

fig = Figure(resolution=PLOT_RES)
ax_k = Axis(fig[1,1], 
    title="k(t), d=$(d)",
    xlabel="t",
    ylabel="k",
    xminorticksvisible = true, 
	xminorgridvisible = true, 
	yminorticksvisible = true, 
	yminorgridvisible = true, 
    xminorticks = IntervalsBetween(10),
	yminorticks = IntervalsBetween(10)
)
ax_k_n = Axis(fig[2,1], 
    title="k(n)",
    xlabel="n",
    ylabel="k",
    xminorticksvisible = true, 
	xminorgridvisible = true, 
	yminorticksvisible = true, 
	yminorgridvisible = true, 
    xminorticks = IntervalsBetween(10),
	yminorticks = IntervalsBetween(10)
)

for i in 1:N
    lines!(ax_k, t_sol, k_sol[i,:], label="$(i)")
end
scatter!(ax_k_n, 1:N, k_final)

if N < 10
    axislegend(ax_k, position=:rt)
end

save(savingpath, fig, px_per_unit=PLOT_PX_PER_UNIT_PNG)