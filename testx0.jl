using Flight

struct EmptyComp <: AbstractAirframeComponent end

function HybridSystem(cmp::EmptyComp, ẋ = get_x0(cmp), x = get_x0(cmp),
                        y = get_y0(cmp), u = get_u0(cmp), d = get_d0(cmp),
                        t = Ref(0.0))
    params = cmp #params is the component itself
    subsystems = nothing #no subsystems to define
    HybridSystem{map(typeof, (cmp, x, y, u, d, params, subsystems))...}(ẋ, x, y, u, d, t, params, subsystems)
end
