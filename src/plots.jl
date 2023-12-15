
using CairoMakie

function beautiful_axis(fit_part; title=nothing, xlabel=nothing, ylabel=nothing, is_x_log=false, is_y_log=false)
    ax = Axis(fit_part, 
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
        xscale=is_x_log ? log10 : identity,
        yscale=is_y_log ? log10 : identity,
        xminorticksvisible = true, 
        xminorgridvisible = true, 
        yminorticksvisible = true, 
        yminorgridvisible = true, 
        xminorticks = IntervalsBetween(10),
        yminorticks = IntervalsBetween(10)
    )
    return ax
end