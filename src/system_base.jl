mutable struct FitzHugh_Nagumo_base_system <: FitzHugh_Nagumo_system 
    p # Постоянные параметры системы
    U₀ # Начальные условия системы

    alg # Алгоритм интегрирования системы
    reltol # Относительная точность интегрирования системы
    abstol # Абсолютная точность интегрирования системы
    maxiters # Максимальное количество итераций
    check_success # Проверка  
    
    sol # Итог численного интегрирования

    FitzHugh_Nagumo_base_system(p, U₀) = begin
        new(p, U₀, ALG, RELTOL, ABSTOL, MAXITERS, false, NaN)
    end
end

""" u = [x, y]\\
p = [a, b, ε] """
function FitzHugh_Nagumo_base_model(du, u, p, t)
    x, y = u
    a, b, ε = p

    du[1] = -x*(x-a)*(x-b) - y
    du[2] = ε*(x-a)
end

""" U₀ = [x, y]\\
param = [a, b, ε] """
function system_integrate!(sys::FitzHugh_Nagumo_base_system, t_span)
    f = FitzHugh_Nagumo_base_model
    prob = ODEProblem(f, sys.U₀, t_span, sys.p)
    sys.sol = solve(prob, sys.alg(), reltol=sys.reltol, abstol=sys.abstol, maxiters=sys.maxiters)
end