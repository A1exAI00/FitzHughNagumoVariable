#=
Из параграфа Chaotic and Regular Patterns in Two-Dimentional Lattices of Coupled Bistable Units
была взята модель бистабильного элемента kⱼ.
В ходе создания цепочки из таких элементов они приходят к разным состояниям равновесия, 
и в зависимости от величины их связи d можно добиться сильного "расщипления". 

В данном скрипте численно интегрируется N осцилляторов ФитцХью-Нагумо, 
каждый из которых после времени t₀ подстраивает свою переменную bⱼ под величину kⱼ - 
величину бистабильного элемента. 
Показываются 4 графика:
1. Эпюры колебаний xⱼ(t)
2. Эпюры bⱼ(t)
3. Эпюры kⱼ(t)
4. Эпюры Tⱼ(t) - период колебаний xⱼ(t)
=#

#########################################################################################

using DrWatson
@quickactivate "FitzHughNagumoVariable"

include(srcdir("plots.jl"))

include(srcdir("FitzHughNagumoVariable.jl"))
using .FitzHughNagumoVariable

#########################################################################################

# Настройки генерируемого графика
PLOT_RES = (1000, 1600)
PLOT_FILENAME = "05-system-2-$(time_ns()).png"
PLOT_SAVE_PATH = plotsdir(PLOT_FILENAME)
PLOT_PX_PER_UNIT_PNG = 2

#########################################################################################

# Постоянные параметры системы
N, ε, μ, d = 2, 1e-3, 1e-1, 1e-2
step_time_scale = 100
a = 0.5
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
U₀ = [x₀ⱼ..., y₀ⱼ..., b₀ⱼ..., k₀ⱼ...]
t_span = [t_min, t_max]

println("system_param=$(system_param)")
println("U₀=$(U₀)")

#########################################################################################

sys = FitzHugh_Nagumo_2_system(system_param, U₀)
system_integrate!(sys, t_span)
sol = sys.sol
t_sol = sol.t

x_sol = Any[]
b_sol = Any[]
k_sol = Any[]
T = Any[]
for i in 1:N
    push!(x_sol, sol[i,:])
    push!(b_sol, sol[2N+i,:])
    push!(k_sol, sol[3N+i,:])
    push!(T, measure_period(sol[i,:], t_sol))
end

k_final = [k_sol[i][end] for i in 1:N]
T_final = [T[i][end] for i in 1:N]

println("k_final = $(k_final)")
println("T_final = $(T_final)")

#########################################################################################

fig = Figure(size=PLOT_RES)
ax_x = beautiful_axis(fig[1,1], title="x(t)", xlabel="t", ylabel="x")
ax_b = beautiful_axis(fig[2,1], title="b(t)", xlabel="t", ylabel="b")
ax_k = beautiful_axis(fig[3,1], title="k(t)", xlabel="t", ylabel="k")
ax_T = beautiful_axis(fig[4,1], title="T(t)", xlabel="t", ylabel="T")

hlines!.((ax_x, ax_b, ax_k, ax_T), 0.0, color=:black)
vlines!.((ax_x, ax_b, ax_k, ax_T), 0.0, color=:black)

for i in 1:N
    lines!(ax_x, t_sol, x_sol[i], label="$(i)")
    lines!(ax_b, t_sol, b_sol[i], label="")
    lines!(ax_k, t_sol, k_sol[i], label="")
    lines!(ax_T, t_sol, T[i], label="T=$(T_final[i])")
end

(N ≤ 10) && axislegend(ax_x, position=:rt)

save(PLOT_SAVE_PATH, fig, px_per_unit=PLOT_PX_PER_UNIT_PNG)