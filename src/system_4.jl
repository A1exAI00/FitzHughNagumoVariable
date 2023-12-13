mutable struct FitzHugh_Nagumo_4_system <: FitzHugh_Nagumo_system 
    p # Постоянные параметры системы
    U₀ # Начальные условия системы

    alg # Алгоритм интегрирования системы
    reltol # Относительная точность интегрирования системы
    abstol # Абсолютная точность интегрирования системы
    maxiters # Максимальное количество итераций
    check_success # Проверка
    is_sol_needed # Надо ли численно интегрировать всю систему
    
    k_sol # Итог численного интегрирования цепочки kⱼ
    sols # Итог численного интегрирования

    FitzHugh_Nagumo_4_system(p, U₀) = begin 
        new(p, U₀, ALG, RELTOL, ABSTOL, MAXITERS, false, true, NaN, NaN)
    end
    FitzHugh_Nagumo_4_system(p, U₀, is_sol_needed) = begin 
        new(p, U₀, ALG, RELTOL, ABSTOL, MAXITERS, false, is_sol_needed, NaN, NaN)
    end
end

""" u = [x, y, b] \\
p = [N, t₀, Δt, ε, μ, a, k] """
function FitzHugh_Nagumo_4_single(du,u,p,t)
    N, t₀, Δt, ε, μ, a, k = p
    x, y, b = u

    du[1] = -x*(x-a)*(x-b) - y
    du[2] = ε*(x-a)
    du[3] = 1/μ * step_func_smooth(t, t₀, Δt) * (k-b)

    nothing
end

""" u = [kⱼ...] \\
p = [N, d, f] """
function FitzHugh_Nagumo_4_k(du,u,p,t)
    N, d, f = p[1:3]
    kⱼ = u[1:end]

    for j in 1:N
        if j == 1
            du[j] = f(kⱼ[j]) + d*(kⱼ[j] - 2*kⱼ[j] + kⱼ[j+1])
        elseif j == N
            du[j] = f(kⱼ[j]) + d*(kⱼ[j-1] - 2*kⱼ[j] + kⱼ[j])
        else
            du[j] = f(kⱼ[j]) + d*(kⱼ[j-1] - 2*kⱼ[j] + kⱼ[j+1])
        end
    end
    nothing
end

""" U₀ = [xⱼ..., yⱼ..., bⱼ..., kⱼ...] \\
param = [N, t₀, Δt, ε, μ, d, f, aⱼ...] """
function system_integrate!(sys::FitzHugh_Nagumo_4_system, t_span)
    k_f = FitzHugh_Nagumo_4_k
    N = sys.p[1]
    k₀ⱼ = sys.U₀[3N+1:4N]; k_p = sys.p[[1,6,7]]

    k_prob = ODEProblem(k_f, k₀ⱼ, t_span, k_p)
    sys.k_sol = solve(k_prob, CVODE_BDF(), reltol=sys.reltol, abstol=sys.abstol, maxiters=sys.maxiters)
    k_final = sys.k_sol[:,end]

    if !(sys.is_sol_needed)
        return nothing
    end

    f = FitzHugh_Nagumo_4_single

    sols = []
    for i in 1:N
        U₀ = sys.U₀[[0N+i, 1N+i, 2N+i]]
        param = (sys.p[[1:5..., 7+i]]..., k_final[i])

        if i == 1 
            prob_precompile = ODEProblem(f, U₀, [0.0, 0.1], param)
            solve(prob_precompile, sys.alg(), reltol=sys.reltol, abstol=sys.abstol, maxiters=sys.maxiters)
        end

        prob = ODEProblem(f, U₀, t_span, param)
        sol = solve(prob, sys.alg(), reltol=sys.reltol, abstol=sys.abstol, maxiters=sys.maxiters)
        push!(sols, sol)
    end
    sys.sols = sols
end