mutable struct FitzHugh_Nagumo_2_system <: FitzHugh_Nagumo_system 
    p # Постоянные параметры системы
    U₀ # Начальные условия системы

    alg # Алгоритм интегрирования системы
    reltol # Относительная точность интегрирования системы
    abstol # Абсолютная точность интегрирования системы
    maxiters # Максимальное количество итераций
    check_success # Проверка  
    
    sol # Итог численного интегрирования

    FitzHugh_Nagumo_2_system(p, U₀) = begin 
        new(p, U₀, CVODE_BDF, 1e-10, 1e-10, MAXITERS, false, NaN)
    end
end

""" u = [xⱼ..., yⱼ..., bⱼ..., kⱼ...] \\
p = [N, t₀, Δt, ε, μ, d, f, aⱼ...] """
function FitzHugh_Nagumo_2(du,u,p,t)
    N, t₀, Δt, ε, μ, d, f = p[1:7]
    aⱼ = p[8:end]

    xⱼ = u[0N+1:1N]
    yⱼ = u[1N+1:2N]
    bⱼ = u[2N+1:3N]
    kⱼ = u[3N+1:4N]

    for j in 1:N
        x, y, b, a = xⱼ[j], yⱼ[j], bⱼ[j], aⱼ[j]

        du[j] = -x*(x-a)*(x-b) - y
        du[N+j] = ε*(x-a)
        du[2N+j] = 1/μ * step_func_smooth(t, t₀, Δt) * (kⱼ[j]-b)

        if j == 1
            du[3N+j] = f(kⱼ[j]) + d*(kⱼ[j] - 2*kⱼ[j] + kⱼ[j+1])
        elseif j == N
            du[3N+j] = f(kⱼ[j]) + d*(kⱼ[j-1] - 2*kⱼ[j] + kⱼ[j])
        else
            du[3N+j] = f(kⱼ[j]) + d*(kⱼ[j-1] - 2*kⱼ[j] + kⱼ[j+1])
        end
    end
    nothing
end

""" u = [xⱼ..., yⱼ..., bⱼ..., kⱼ...] \\
p = [N, t₀, Δt, ε, μ, d, f, aⱼ...] """
function FitzHugh_Nagumo_2_jac(J,u,p,t)
    N, t₀, Δt, ε, μ, d, f = p[1:7]
    aⱼ = p[8:end]

    xⱼ = u[0N+1:1N]
    yⱼ = u[1N+1:2N]
    bⱼ = u[2N+1:3N]
    kⱼ = u[3N+1:4N]
    
    for i in 1:N
        x, y, b, a, k = xⱼ[i], yⱼ[i], bⱼ[i], aⱼ[i], kⱼ[i]

        # dot x 
        J[i,i] = -3x^2 + 2x*(a+b) - a*b
        J[i,N+i] = -1
        J[i,2N+i] = x*(x-a)

        # dot y
        J[N+i,i] = ε

        # dot b 
        J[2N+i,2N+i] = -1/μ * step_func_smooth(t, t₀, Δt)
        J[2N+i,3N+i] = 1/μ * step_func_smooth(t, t₀, Δt)

        # dot k
        if (i == 1) || (i == N)
            J[3N+i,3N+i] = (-3k^2 + 9k - 6.5) - d
        else
            J[3N+i,3N+i] = (-3k^2 + 9k - 6.5) - 2d
        end
    end

    for i in 1:N-1
        J[3N+i,3N+i+1] = d
        J[3N+i+1,3N+i] = d
    end
    nothing
end

""" U₀ = [xⱼ..., yⱼ..., bⱼ..., kⱼ...] \\
param = [N, t₀, Δt, ε, μ, d, f, aⱼ...] """
function system_integrate!(sys::FitzHugh_Nagumo_2_system, t_span)
    f = FitzHugh_Nagumo_2
    f_jac = FitzHugh_Nagumo_2_jac
    prob = ODEProblem(ODEFunction(f; jac=f_jac), sys.U₀, t_span, sys.p)
    sys.sol = solve(prob, sys.alg(), reltol=sys.reltol, abstol=sys.abstol, maxiters=sys.maxiters)
end