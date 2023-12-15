#=
Из параграфа Chaotic and Regular Patterns in Two-Dimentional Lattices of Coupled Bistable Units
была взята модель бистабильного элемента kⱼ.

В данном скрипте интегрируется цепочка из N бистабильных элементов. Каждый элемент в начальный 
момент времени имеет случайное значение.
=#

#########################################################################################

using DrWatson
@quickactivate "FitzHughNagumoVariable"
using Distributions: Uniform

include(srcdir("plots.jl"))
include(srcdir("FitzHughNagumoVariable.jl"))
using .FitzHughNagumoVariable

#########################################################################################

# Настройки генерируемого графика
# TODO refactor naming and saving names
PLOT_RES = (600, 500)
PLOT_SAVING_DIR = plotsdir(); println(pwd())
PLOT_FILENAME = "09-k_chain-VI-"
savingpath = joinpath(PLOT_SAVING_DIR, PLOT_FILENAME*"$(time_ns()).png")
PLOT_PX_PER_UNIT_PNG = 2

#########################################################################################

# Постоянные параметры системы
N = 100
ε = 1e-3
μ = 1e-1
a = 0.5
d = 0.07
step_time_scale = 100
c = 0.65

f(x,c) = x*(x-1)*(c-x)
f(x) = f(x,c)

# Начальные условия системы
x₀, y₀, b₀ = 0.0, 0.0, 1.5
k₀_min, k₀_max = 0.0, 1.0

# Время интегрирования
t_min = 0.0
t_max = 100.0
t₀ = (t_max-t_min)/10

#########################################################################################

aⱼ = a .* ones(N)
system_param = (N, t₀, step_time_scale, ε, μ, d, f, aⱼ...)
x₀ⱼ = x₀ .* ones(N)
y₀ⱼ = y₀ .* ones(N)
b₀ⱼ = b₀ .* ones(N)
k₀ⱼ = rand(Uniform(k₀_min, k₀_max), N)
U₀ = [x₀ⱼ..., y₀ⱼ..., b₀ⱼ...,k₀ⱼ...]
t_span = (t_min, t_max)

#########################################################################################

sys = FitzHugh_Nagumo_3_system(system_param, U₀, false)
system_integrate!(sys, t_span)
k_sol = sys.k_sol
t_sol = k_sol.t

k_final = k_sol[:,end]

#########################################################################################

fig = Figure(size=PLOT_RES)
ax_k = beautiful_axis(fig[1,1], title="k(t), d=$(d), c=$(c), kⱼ₀∈($(k₀_min), $(k₀_max))", xlabel="t", ylabel="k")
ax_k_n = beautiful_axis(fig[2,1], title="k(n)", xlabel="n", ylabel="k")

for i in 1:N
    lines!(ax_k, t_sol, k_sol[i,:], label="$(i)")
end
scatter!(ax_k_n, 1:N, k_final)

# TODO refactor
if N < 10
    axislegend(ax_k, position=:rt)
end

save(savingpath, fig, px_per_unit=PLOT_PX_PER_UNIT_PNG)