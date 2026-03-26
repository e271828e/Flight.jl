module FlightPhysicsControlSystemsBaseExt

using ControlSystemsBase
using FlightPhysics.Linearization

function ControlSystemsBase.ss(lss::LinearizedSS)
    ControlSystemsBase.ss(lss.A, lss.B, lss.C, lss.D)
end

function Linearization.LinearizedSS(mdl::ControlSystemsBase.StateSpace{ControlSystemsBase.Continuous, <:AbstractFloat})
    (; A, B, C, D, nx, nu, ny) = mdl
    ẋ0 = zeros(nx); x0 = zeros(nx); u0 = zeros(nu); y0 = zeros(ny)
    LinearizedSS(; ẋ0, x0, u0, y0, A, B, C, D)
end


end
