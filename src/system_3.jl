mutable struct FitzHugh_Nagumo_3_system <: FitzHugh_Nagumo_system 
    p # Постоянные параметры системы
    U₀ # Начальные условия системы

    alg # Алгоритм интегрирования системы
    reltol # Относительная точность интегрирования системы
    abstol # Абсолютная точность интегрирования системы
    maxiters # Максимальное количество итераций
    check_success # Проверка
    is_sol_needed # Надо ли численно интегрировать всю систему
    
    k_sol # Итог численного интегрирования цепочки kⱼ
    sol # Итог численного интегрирования

    FitzHugh_Nagumo_3_system(p, U₀) = begin 
        new(p, U₀, ALG, RELTOL, ABSTOL, MAXITERS, false, true, NaN, NaN)
    end
    FitzHugh_Nagumo_3_system(p, U₀, is_sol_needed) = begin 
        new(p, U₀, ALG, RELTOL, ABSTOL, MAXITERS, false, is_sol_needed, NaN, NaN)
    end
end

""" u = [xⱼ..., yⱼ..., bⱼ...] \\
p = [N, t₀, Δt, ε, μ, aⱼ..., kⱼ...] """
function FitzHugh_Nagumo_3(du,u,p,t)
    N, t₀, Δt, ε, μ = p[1:5]
    aⱼ = p[0N+6:1N+5]
    kⱼ = p[1N+6:end]

    xⱼ = u[0N+1:1N]
    yⱼ = u[1N+1:2N]
    bⱼ = u[2N+1:3N]

    for j in 1:N
        x, y, b, a, k = xⱼ[j], yⱼ[j], bⱼ[j], aⱼ[j], kⱼ[j]

        du[j] = -x*(x-a)*(x-b) - y
        du[N+j] = ε*(x-a)
        du[2N+j] = 1/μ * step_func_smooth(t, t₀, Δt) * (k-b)
    end
    nothing
end

""" u = [kⱼ...] \\
p = [N, d, f] """
function FitzHugh_Nagumo_3_k(du,u,p,t)
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
function system_integrate!(sys::FitzHugh_Nagumo_3_system, t_span)
    k_f = FitzHugh_Nagumo_3_k
    N = sys.p[1]
    k₀ⱼ = sys.U₀[3N+1:4N]; k_p = sys.p[[1,6,7]]

    k_prob = ODEProblem(k_f, k₀ⱼ, t_span, k_p)
    sys.k_sol = solve(k_prob, CVODE_BDF(), reltol=sys.reltol, abstol=sys.abstol, maxiters=sys.maxiters)
    k_final = sys.k_sol[:,end]

    if !(sys.is_sol_needed)
        return nothing
    end

    f = FitzHugh_Nagumo_3
    U₀ = sys.U₀[1:3N]
    param = (sys.p[[1:5..., 8:N+7...]]..., k_final...)
    
    prob_precompile = ODEProblem(f, U₀, [0.0, 0.1], Tuple(param))
    solve(prob_precompile, sys.alg(), reltol=sys.reltol, abstol=sys.abstol, maxiters=sys.maxiters)

    prob = ODEProblem(f, U₀, t_span, Tuple(param))
    sys.sol = solve(prob, sys.alg(), reltol=sys.reltol, abstol=sys.abstol, maxiters=sys.maxiters)
end