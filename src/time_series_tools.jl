#=
A collection of function, that manipulate time series.
=#

########################################################################

""" Гладкая альтернатива функции Хевисайда с характерным временным масштабом Δt """
function step_func_smooth(t, t_step, Δt)
    k = 1/Δt * log(3)
    0.5 * (tanh(k*(t-t_step)) + 1)
end

########################################################################

"""
`tail_fraction` : отношение хвоста `seq` ко всему `seq` \\
`0≤tail_fraction≤1`"""
function mean_of_tail(seq, tail_fraction)
    tail_start_index = round(Int, tail_fraction*length(seq))
    tail_seq = seq[tail_start_index:end]
    return mean(tail_seq)
end

""" bounds=(min,max) """
function is_out_of_bounds(seq, bounds)
    min, max = bounds
    return any(i->(i<min || i>max), seq)
end

""" bounds=(min,max) """
function is_in_bounds(seq, bounds)
    return !(is_out_of_bounds(seq, bounds))
end

########################################################################

function indexes_of_maxes(seq)
    indexes = []
    for i in 2:length(seq)-1
        x₁, x₂, x₃ = seq[i-1:i+1]
        if (x₁<x₂ && x₃<x₂) push!(indexes, i) end
    end
    return indexes
end

function indexes_of_mins(seq)
    indexes = []
    for i in 2:length(seq)-1
        x₁, x₂, x₃ = seq[i-1:i+1]
        if (x₁>x₂ && x₃>x₂) push!(indexes, i) end
    end
    return indexes
end

function times_of_maxes(seq, t_seq)
    indexes = indexes_of_maxes(seq)
    return t_seq[indexes]
end

function times_of_mins(seq, t_seq)
    indexes = indexes_of_mins(seq)
    return t_seq[indexes]
end

"""Measures period of `seq`

# Parameters

`seq` : array

`t_seq` : array of time

# Algorithm

Function calculates mean time difference between maximums, 
mean time difference between minimums and avereges them.
Throws error if there is less then 2 of maximums or minimums."""
function measure_period_avg(seq, t_seq)
    times_max = times_of_maxes(seq, t_seq)
    times_min = times_of_mins(seq, t_seq)
    (length(times_max) < 2) && (length(times_min) < 2) && return NaN
    T = Vector{Float64}()
    (length(times_max) < 2) || push!(T, mean(diff(times_max))) 
    (length(times_min) < 2) || push!(T, mean(diff(times_min))) 
    return mean(T)
end

function measure_period(seq, t_seq)
    times_max = times_of_maxes(seq, t_seq)
    N_maxes = length(times_max)

    period = zeros(length(seq))

    for i in eachindex(seq)
        t = t_seq[i]
        if t ≤ times_max[1]
            period[i] = times_max[2] - times_max[1]
            continue
        elseif times_max[end] ≤ t 
            period[i] = times_max[end] - times_max[end-1]
            continue
        end
        for j in 1:N_maxes-1
            if times_max[j] ≤ t ≤ times_max[j+1]
                period[i] = times_max[j+1] - times_max[j]
                break
            end
        end
    end
    return period
end

# function generate_square_itp(itp_param)
#     t₀, t₁, tₙ, T, D, A, ξ = itp_param
#     t_range = range(t₀, t₁, tₙ)
#     A₁, A₂ = 2*A*(1-D), -2*A*D
#     func(t) = (T*(1-D)/2 ≤ t ≤ T*(1+D)/2) ? A₁ : A₂
#     func_noisy(t) = func(t) + ξ*(2*rand()-1)
#     signal_noisy = func_noisy.(t_range)
#     return linear_interpolation(t_range, signal_noisy)
# end

# function generate_cos_itp(itp_param)
#     t₀, t₁, tₙ, N, ξ = itp_param[1:5]
#     Ω_teach = itp_param[0N+6:1N+5]
#     A_teach = itp_param[1N+6:2N+5]
#     Φ_teach = itp_param[2N+6:3N+5]
#     t_range = range(t₀, t₁, tₙ)
#     func(t) = sum(A_teach.*cos.(Ω_teach.*t .+ Φ_teach))
#     func_noisy(t) = func(t) + ξ*(2*rand()-1)
#     signal_noisy = func_noisy.(t_range)
#     return linear_interpolation(t_range, signal_noisy)
# end