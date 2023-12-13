#=
Базовые функции, необходимые для интегрирования модели 
ФитцХью-Нагумо и адаптивной модели ФитцХью-Нагумо, исследуемой в данной работе
=#

#########################################################################################



#########################################################################################

# Функция медленного движения
f_slow(x,a,b) = -x*(x-a)*(x-b)

# X координаты точек A,B,C,D
xA(a,b) = 1/3 * (a+b + sqrt(a^2+b^2-a*b))
xB(a,b) = 1/3 * (a+b-2*sqrt(a^2+b^2-a*b))
xC(a,b) = 1/3 * (a+b - sqrt(a^2+b^2-a*b))
xD(a,b) = 1/3 * (a+b+2*sqrt(a^2+b^2-a*b))

# Y координаты точек A,B,C,D
yA(a,b) = 1/27 * (2*(a^3+b^3) - 3*a*b*(a+b) + 2*(a^2+b^2-a*b)^(3/2))
yB(a,b) = yA(a,b)
yC(a,b) = 1/27 * (2*(a^3+b^3) - 3*a*b*(a+b) - 2*(a^2+b^2-a*b)^(3/2))
yD(a,b) = yC(a,b)

# Подынтегральное выражение для периода
G_int_2(x,a,b,ε) = 1/ε/(x-a) * (-3x^2+2x*(a+b)-a*b) # Сходится с интегралом с G_int_2(x,a,b,ε)
G(x,p) = G_int_2(x,p...)

# Интегрирование периода 
# (Индексация нужна из-за того, что метод численного интегрирования 
# возвращает 0-мерный массив, и чтобы получить из него float, берем первый (и единственный) 
# элемент этого массива)
T_BC_prob(a,b,ε) = IntegralProblem(G, xB(a,b), xC(a,b), (a,b,ε))
T_DA_prob(a,b,ε) = IntegralProblem(G, xD(a,b), xA(a,b), (a,b,ε))
T_BC_solve(a,b,ε) = abs(solve(T_BC_prob(a,b,ε), QuadGKJL())[1])
T_DA_solve(a,b,ε) = abs(solve(T_DA_prob(a,b,ε), QuadGKJL())[1])
T_Σ_int(a,b,ε) = T_BC_solve(a,b,ε) + T_DA_solve(a,b,ε)

# Аналитически решенный интеграл с G_int_2(x,a,b,ε)
T_analitic_(a,b,ε,x1,x2) = 1/ε * (
    -3/2*(x2^2-x1^2) +
    (x2-x1)*(2b-a) + 
    log(abs((x2-a)/(x1-a)))*(a*b-a^2)
)
T_BC_analitic(a,b,ε) = T_analitic_(a,b,ε, xB(a,b),xC(a,b))
T_DA_analitic(a,b,ε) = T_analitic_(a,b,ε, xD(a,b),xA(a,b))
T_Σ_analitic(a,b,ε) = T_BC_analitic(a,b,ε) + T_DA_analitic(a,b,ε)

#########################################################################################

abstract type FitzHugh_Nagumo_system end

include("system_base.jl")
include("system_1.jl")
include("system_2.jl")
include("system_3.jl")
include("system_4.jl")