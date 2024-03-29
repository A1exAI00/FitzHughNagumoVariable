#=
В данном скрипте реализовано то же самое, что 05-model2_learning.jl, 
но разделяя систему на 2: осцилляторы и эволюция k во времени.
=#

#########################################################################################

using DrWatson
@quickactivate "FitzHughNagumoVariable"

include(srcdir("plots.jl"))

include(srcdir("FitzHughNagumoVariable.jl"))
using .FitzHughNagumoVariable

#########################################################################################

# Настройки генерируемого графика
PLOT_RES = (1000, 1000)
PLOT_FILENAME = "07-system-3-$(time_ns()).png"
PLOT_SAVE_PATH = plotsdir(PLOT_FILENAME)
PLOT_PX_PER_UNIT_PNG = 2

#########################################################################################

# Постоянные параметры системы
N = 50 # Количество осцилляторов
a = 0.5
ε = 1e-3
μ = 1e-1
d = 1e-1 # Коэффициент связи между элементами цепочки
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

sys = FitzHugh_Nagumo_3_system(system_param, U₀)
system_integrate!(sys, t_span)
sol = sys.sol
t_sol = sol.t

x_sol, b_sol = Any[], Any[]
for i in 1:N
    push!(x_sol, sol[i,:])
    push!(b_sol, sol[2N+i,:])
end

k_sol = sys.k_sol
k_t_sol = k_sol.t

#########################################################################################

fig = Figure(size=PLOT_RES)
ax_x = beautiful_axis(fig[1,1], title="x(t)", xlabel="t", ylabel="x")
ax_b = beautiful_axis(fig[2,1], title="b(t)", xlabel="t", ylabel="b")
ax_k = beautiful_axis(fig[3,1], title="k(t), d=$(d)", xlabel="t", ylabel="k")
ax_b_n = beautiful_axis(fig[4,1], title="b(n), d=$(d)", xlabel="n", ylabel="b")

# axis
hlines!.((ax_x, ax_b, ax_k, ax_b_n), 0.0, color=:black)
vlines!.((ax_x, ax_b, ax_k, ax_b_n), 0.0, color=:black)

for i in 1:N
    lines!(ax_x, t_sol, x_sol[i], label="$(i)")
    lines!(ax_b, t_sol, b_sol[i], label="")
    lines!(ax_k, k_t_sol, k_sol[i,:], label="")
end
scatter!(ax_b_n, 1:N, k_sol[:,end])

(N < 10) && axislegend(ax_x, position=:rt)

save(PLOT_SAVE_PATH, fig, px_per_unit=PLOT_PX_PER_UNIT_PNG)