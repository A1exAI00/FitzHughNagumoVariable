module FitzHughNagumoVariable

    using Statistics: mean
    using OrdinaryDiffEq, StaticArrays, Sundials
    using Integrals
    using Dates

    include("misc_tools.jl")
    include("time_series_tools.jl")
    include("oscillator_models.jl")

    #########################################################################################

    # const ALG = CVODE_BDF
    const ALG = Rodas5P
    const RELTOL::Float64, ABSTOL::Float64 = 1e-5, 1e-5
    const MAXITERS::Int64 = Int(1e7)

    #########################################################################################

    # из oscillator_models.jl
    export f_slow, xA, xB, xC, xD, yA, yB, yC, yD, T_Σ_analitic
    export step_func_smooth
    export system_integrate!
    export FitzHugh_Nagumo_base_system
    export FitzHugh_Nagumo_1_system
    export FitzHugh_Nagumo_2_system
    export FitzHugh_Nagumo_3_system
    export FitzHugh_Nagumo_4_system

    # из time_series_tools.jl
    export mean_of_tail, is_out_of_bounds, is_in_bounds
    export indexes_of_maxes, indexes_of_mins
    export times_of_maxes, times_of_mins
    export measure_period_avg, measure_period

    # из misc_tools.jl
    export elapsed_time_string, get_easy_time

end