β = ISA_layers[1].β

compute_δ_p(p) = (p/p_std) ^ (1+β*R/(2*g_std)) #computes normalized δ

function compute_δ_h(h) #computes normalized δ
    @unpack p, T = ISAData(AltGeop(h))
    p / p_std / √(T / T_std)
end

@assert compute_δ_h(1345) ≈ compute_δ_p(ISAData(AltGeop(1345)).p)

#our inputs are p, T

function f_cont!(n, MAP, δ)

    #ω es propiedad del motor? no, porque no puede calcularla solo con sus
    #datos. para eso necesita el momento de inercia, que no depende solo de el.
    #no, el motor no tiene ω como estado. la cuestion es: el motor deberia ser
    #System o ni siquiera? quiza si, porque lo que si tiene como estado es
    #on/off. a lo mejor debe tener un estado discreto

    #el motor debe ser System

    #engine inputs: throttle
    #engine outputs: throttle, generated_torque, power, fuel consumption, throttle, MAP

    #aunque realmente si podriamos hacer que omega fuera del motor. si le
    #pasamos el torque y el momento de inercia equivalente que esta moviendo. y
    #eso ya si lo calcula su parent System powerplant. llamando a la helice.
    #entonces la ecuacion de momento cinetico axial la puede aplicar el propio
    #motor. tiene esto sentido??

    #en este caso a f_cont! habria que pasarle solo el torque y el momento de
    #inercia en el eje. que se encarga el parent System
    #f_cont!(sys::System{<:PistonEngine}, air::AirData, T_load::Real, Ixx_load::Real)
    #T_load deberia ser negativa en x. O sea, la gearbox se las tiene que apanar
    #para que el motor siempre vea un par resistente CCW, porque un motor es
    #siempre CW.

    #get x, from x omega, omega to n, from n and n_rated, n̄.
    #get air data. from p, delta. from n̄ and delta, M̃_wot

    #compute the wide-open throttle MAP for the given RPMs and altitude
    @show Mbar_wot = interp_Mbar_wot(nbar, δ)
    idle_MAP_ratio = 0.4
    #this can be tuned so that the engine idles at appropriate RPMs with the
    #chosen propeller
    @show Mbar = Mbar_wot * (idle_MAP_ratio + thr * (1 - idle_MAP_ratio))

    #δ at which our Mbar would be Mbar_wot

    Pbar_ISA_std = interp_Pbar_ISA_std(nbar, Mbar)

    @show δ_wot = interp_δ_wot(nbar, Mbar)
    Pbar_ISA_wot = interp_Pbar_ISA_wot(nbar, δ_wot)
    @show P_ISA_std = Pbar_ISA_std * P_rated
    @show P_ISA_wot = Pbar_ISA_wot * P_rated

    #when p_wot is close to p_std with MAP = MAP_wot (thr = 1),
    #p_wot(MAP_wot(p_std)) = p_std and P̃_wot = P̃_std. we need to avoid the
    #division by zero
    @show abs(δ_wot - 1)
    if abs(δ_wot - 1) < 1e-3
        @show "Hi"
        Pbar_ISA = Pbar_ISA_std
    else
        Pbar_ISA = Pbar_ISA_std + (Pbar_ISA_wot - Pbar_ISA_std) / (δ_wot - 1) * (δ - 1)
        @show Pbar_ISA_std
    end

    @show P_ISA = Pbar_ISA * P_rated
    @show MAP = Mbar * MAP_rated
    @show δ = Mbar * MAP_rated

    P =

    return max(0, P_ISA)

    #we need to return manifold pressure

end