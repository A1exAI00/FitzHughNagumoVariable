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
using Statistics: mean, std

include(srcdir("plots.jl"))
include(srcdir("FitzHughNagumoVariable.jl"))
using .FitzHughNagumoVariable

#########################################################################################

# Настройки генерируемого графика
# TODO refactor naming and saving names
PLOT_RES = (500, 500)
PLOT_SAVING_DIR = plotsdir(); println(pwd())
PLOT_FILENAME = "10-Dch_region-VI-"
savingpath = joinpath(PLOT_SAVING_DIR, PLOT_FILENAME*"$(time_ns()).png")
PLOT_PX_PER_UNIT_PNG = 2

#########################################################################################

# Постоянные параметры системы
N = 40
ε = 1e-3
μ = 1e-1
a = 0.5
step_time_scale = 100
d_min, d_max, N_d = 0.0, 0.15, 100
c_min, c_max, N_c = 0.0, 1.0, 100
f(x,c) = x*(x-1)*(c-x)

N_sims = 30

# Начальные условия системы
x₀, y₀, b₀ = 0.0, 0.0, 1.5
k₀_min, k₀_max = -1.0, 1.0
# k₀_min, k₀_max = 0.0, 1.0

# Время интегрирования
t_min = 0.0
t_max = 10000.0
t₀ = (t_max-t_min)/10

Dch_metric = maximum
# Dch_metric = std

#########################################################################################

aⱼ = a .* ones(N)
c_range = range(c_min, c_max, N_c)
d_range = range(d_min, d_max, N_d)

x₀ⱼ = x₀ .* ones(N)
y₀ⱼ = y₀ .* ones(N)
b₀ⱼ = b₀ .* ones(N)
k₀ⱼ = rand(Uniform(k₀_min, k₀_max), N)
U₀ = [x₀ⱼ..., y₀ⱼ..., b₀ⱼ...,k₀ⱼ...]

t_span = (t_min, t_max)

#########################################################################################

diffs_k = zeros(N_c, N_d)
for i in 1:N_c, j in 1:N_d, k in 1:N_sims
    (rem(k,10)==0) && (rem(j,10)==0) && @show (i,j,k)
    c = c_range[i]
    d = d_range[j]
    f_(x) = f(x,c)
    system_param = (N, t₀, step_time_scale, ε, μ, d, f_, aⱼ...)

    sys = FitzHugh_Nagumo_3_system(system_param, U₀, false)
    system_integrate!(sys, t_span)
    k_sol = sys.k_sol
    t_sol = k_sol.t

    k_final = k_sol[:,end]
    diff_k = diff(k_final)
    max_diff_k = Dch_metric(diff_k)
    diffs_k[i,j] += max_diff_k
end

#########################################################################################

fig = Figure(size=PLOT_RES)

# TODO refactor title, line too big
ax_Dch = beautiful_axis(fig[1,1], title="Dch; kⱼ₀∈($(k₀_min), $(k₀_max)), Dch_metric=$(Dch_metric)", xlabel="c", ylabel="d")
xlims!(ax_Dch, c_min, c_max)

contour!(ax_Dch, c_range, d_range, diffs_k)

save(savingpath, fig, px_per_unit=PLOT_PX_PER_UNIT_PNG)