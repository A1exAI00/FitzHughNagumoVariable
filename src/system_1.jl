mutable struct FitzHugh_Nagumo_1_system <: FitzHugh_Nagumo_system 
    p # Постоянные параметры системы
    U₀ # Начальные условия системы

    alg # Алгоритм интегрирования системы
    reltol # Относительная точность интегрирования системы
    abstol # Абсолютная точность интегрирования системы
    maxiters # Максимальное количество итераций
    check_success # Проверка  
    
    sol # Итог численного интегрирования

    FitzHugh_Nagumo_1_system(p, U₀) = begin 
        new(p, U₀, ALG, RELTOL, ABSTOL, MAXITERS, false, NaN)
    end
end

""" u = [x, y, b]\\
p = [a, ε, μ, k, t_crit] """
function FitzHugh_Nagumo_1(du, u, p, t)
    x, y, b = u
    a, ε, μ, k, t_crit = p

    du[1] = -x*(x-a)*(x-b) - y
    du[2] = ε*(x-a)
    du[3] = 1/μ * step_func_smooth(t, t_crit, 50)*(k-b)
end

""" U₀ = [x, y, b]\\
param = [a, ε, μ, k, t_crit] """
function system_integrate!(sys::FitzHugh_Nagumo_1_system, t_span)
    f = FitzHugh_Nagumo_1
    prob = ODEProblem(f, sys.U₀, t_span, sys.p)
    sys.sol = solve(prob, sys.alg(), reltol=sys.reltol, abstol=sys.abstol, maxiters=sys.maxiters)
end