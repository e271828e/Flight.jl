using StaticArrays
using StructArrays

function helix(t)
    p = 1
    R = 1
    ω = 1 #only relevant for kinematics

    x = R*cos(ω*t)
    y = R*sin(ω*t)
    z = p*ω*R*t #p = pitch

    return (x = x, y = y, z = z)
end


function plot_helix()

    t_span = 0:0.01:2π |> collect
    data = (helix(t) for t in t_span) |> collect;

    th = TimeHistory(t_span, data)

    # plot(th)
    sa_data = StructArray(data)
    plot(sa_data.x, sa_data.y, sa_data.z)

end